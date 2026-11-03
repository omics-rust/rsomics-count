"""Verify source packaging and the installed count-to-normalization workflow."""

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import tarfile

PRODUCT = 'rsomics-count'
ROOT = Path(__file__).resolve().parents[1]
HOSTS = {'linux-x86_64': 'x86_64-unknown-linux-gnu', 'linux-aarch64': 'aarch64-unknown-linux-gnu',
         'macos-x86_64': 'x86_64-apple-darwin', 'macos-aarch64': 'aarch64-apple-darwin'}


def identity(path):
    body = Path(path).read_bytes()
    return dict(bytes=len(body), sha256=hashlib.sha256(body).hexdigest())


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def package_contract(package):
    parent = Path(package['manifest_path']).parent
    result = {k: v for k, v in package.items() if k not in ('id', 'manifest_path', 'targets')}
    result['targets'] = [{**t, 'src_path': str(Path(t['src_path']).relative_to(parent))} for t in package['targets']]
    return result


def extract_package(source, destination, expected, head, candidate=False):
    destination = Path(destination)
    destination.mkdir(parents=True)
    with tarfile.open(source, 'r:gz') as archive:
        members = archive.getmembers()
        require(len(members) == len({m.name for m in members}), 'duplicate source package member')
        require(len(members) == len({m.name.casefold() for m in members}), 'case-colliding source package member')
        bodies = {}
        for member in members:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and '..' not in path.parts and path.parts[0] == destination.name,
                    'unsafe source package path')
            require(member.isfile() or member.isdir(), 'unsupported source package member')
            if member.isfile():
                name = '/'.join(path.parts[1:])
                mode = expected[name]['mode'] if name in expected else 0o644
                require(name and member.mode == mode, 'source package file mode differs: ' + name)
                bodies[name] = archive.extractfile(member).read()
        require(set(bodies) == set(expected) | {'Cargo.toml.orig', '.cargo_vcs_info.json'},
                'complete source package membership differs')
        require(bodies['Cargo.toml.orig'] == expected['Cargo.toml']['body'], 'original Cargo manifest differs')
        vcs = json.loads(bodies['.cargo_vcs_info.json'])
        require(vcs['git']['sha1'] == head and bool(vcs['git'].get('dirty', False)) == candidate and
                vcs['path_in_vcs'] == '', 'source package VCS identity differs')
        for name, row in expected.items():
            if name != 'Cargo.toml':
                require(bodies[name] == row['body'], 'source package body differs: ' + name)
        for name, body in bodies.items():
            path = destination / name
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open('xb') as stream:
                stream.write(body)
            path.chmod(expected[name]['mode'] if name in expected else 0o644)
    return dict(**identity(source), members=len(bodies), original_members=len(expected),
                original_bodies_modes_match=True, source_commit=head, clean_vcs=not candidate,
                local_uncommitted_candidate=candidate, immutable_registry_archive_equivalent=False)


def source(candidate=False):
    temporary = Path(os.environ['RUNNER_TEMP'])
    delivery = temporary / 'delivery'
    metadata = json.loads((delivery / 'package-metadata.json').read_bytes())
    package, = [p for p in metadata['packages'] if p['name'] == PRODUCT]
    prefix = PRODUCT + '-' + package['version']
    archive = Path(metadata['target_directory']) / 'package' / (prefix + '.crate')
    shutil.copy2(archive, delivery / archive.name)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    require(head == os.environ['GITHUB_SHA'], 'source head differs from actual CI head')
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT) != b''
    require(dirty == candidate, 'clean CI source or explicitly dirty local candidate required')
    expected = {}
    for row in subprocess.check_output(['git', 'ls-tree', '-rz', head], cwd=ROOT).rstrip(b'\0').split(b'\0'):
        header, name = row.split(b'\t')
        mode, kind, blob = header.decode().split()
        name = name.decode()
        require(kind == 'blob' and mode in ('100644', '100755'), 'unexpected source Git type/mode: ' + name)
        body = (ROOT / name).read_bytes()
        require(candidate or hashlib.sha1(b'blob ' + str(len(body)).encode() + b'\0' + body).hexdigest() == blob,
                'source differs from exact Git body: ' + name)
        expected[name] = dict(mode=int(mode, 8) & 0o777, body=body)
    if candidate:
        names = (delivery / 'package-files.txt').read_text().splitlines()
        expected = {name: dict(mode=(ROOT / name).stat().st_mode & 0o777, body=(ROOT / name).read_bytes())
                    for name in names if name not in ('Cargo.toml.orig', '.cargo_vcs_info.json')}
    destination = temporary / 'packaged-source' / prefix
    verification = extract_package(archive, destination, expected, head, candidate)
    unpacked = json.loads(subprocess.check_output(['cargo', 'metadata', '--manifest-path', str(destination / 'Cargo.toml'),
                                                 '--locked', '--offline', '--no-deps', '--format-version', '1']))
    installed, = [p for p in unpacked['packages'] if p['name'] == PRODUCT]
    require(package_contract(installed) == package_contract(package), 'normalized source package metadata differs')
    verification['normalized_metadata_match'] = True
    save(delivery / 'package-verification.json', verification)
    metadata_body = subprocess.check_output(['cargo', 'metadata', '--locked', '--format-version', '1'], cwd=ROOT)
    tree_body = subprocess.check_output(['cargo', 'tree', '--locked', '--offline', '--target', os.environ['EXPECTED_HOST'],
                                        '--edges', 'normal,build', '--prefix', 'none', '--format', '{p}'], cwd=ROOT)
    (delivery / 'native-dependencies.metadata.json').write_bytes(metadata_body)
    (delivery / 'native-dependencies.tree.txt').write_bytes(tree_body)
    feature_body = subprocess.check_output(['cargo', 'tree', '--locked', '--offline', '--target',
                                           os.environ['EXPECTED_HOST'], '--edges', 'normal,build,features'], cwd=ROOT)
    (delivery / 'native-features.tree.txt').write_bytes(feature_body)
    by_name = {(p['name'], p['version']): p for p in json.loads(metadata_body)['packages']}
    enabled = set()
    for line in tree_body.decode().splitlines():
        name, version = re.match(r'([^ ]+) v([^ ]+)', line).groups()
        if name != PRODUCT:
            enabled.add(by_name[name, version]['id'])
    inventory = json.loads((ROOT / 'licenses/provenance.json').read_bytes())
    require(enabled == set(inventory['target_graphs'][os.environ['EXPECTED_HOST']]), 'enabled target dependency inventory differs')
    require(identity(ROOT / 'Cargo.lock')['sha256'] == inventory['lock_sha256'], 'notice inventory lock differs')
    for row in inventory['copied_texts'].values():
        require(identity(ROOT / row['destination']) == {k: row[k] for k in ('bytes', 'sha256')}, 'bundled notice body differs')
    save(delivery / 'native-dependencies.json', dict(target=os.environ['EXPECTED_HOST'], enabled_packages=sorted(enabled),
         inventory_sha256=identity(ROOT / 'licenses/provenance.json')['sha256'],
         metadata_sha256=hashlib.sha256(metadata_body).hexdigest(), tree_sha256=hashlib.sha256(tree_body).hexdigest(),
         feature_tree_sha256=hashlib.sha256(feature_body).hexdigest(),
         enabled_normal_build_graph_verified=True, linked_binary_inventory=False))
    with Path(os.environ['GITHUB_ENV']).open('a') as stream:
        stream.write('PACKAGE_DIRECTORY=' + str(destination) + '\n')
        stream.write('PACKAGE_VERSION=' + package['version'] + '\n')
        stream.write('SOURCE_PACKAGE=' + str(delivery / archive.name) + '\n')


def smoke(binary):
    binary = Path(binary).resolve(strict=True)
    directory = Path(os.environ['RUNNER_TEMP']) / 'delivery/lifecycle'
    directory.mkdir(parents=True)
    shutil.copy2(binary, directory / 'installed-binary')
    fixture = ROOT / 'tests/fixtures'
    shutil.copy2(fixture / 'single.sam', directory / 'library-a.sam')
    (directory / 'library-b.sam').write_bytes(b''.join(line for line in
        (fixture / 'single.sam').read_bytes().splitlines(keepends=True) if not line.startswith(b'rA\t')))
    shutil.copy2(fixture / 'model.gtf', directory / 'genes.gtf')
    for name in ('projection-paired.sam', 'projection-paired.saf'):
        shutil.copy2(fixture / name, directory / name)
    records = []

    def run(arguments, *, stdin=None, status=0):
        command = [str(binary), *map(str, arguments)]
        index = len(records)
        row = dict(argv=command, stdin=None if stdin is None else dict(path=str(stdin), **identity(stdin)))
        records.append(row)
        save(directory / 'commands.json', records)
        with (directory / f'{index}.stdout').open('xb') as out, (directory / f'{index}.stderr').open('xb') as err:
            result = subprocess.run(command, cwd=directory, input=None if stdin is None else Path(stdin).read_bytes(),
                                    stdout=out, stderr=err, env=dict(os.environ, NO_COLOR='1'), timeout=120)
        row['returncode'] = result.returncode
        save(directory / 'commands.json', records)
        require(result.returncode == status, 'installed command exit differs: ' + str(command))
        return directory / f'{index}.stdout', directory / f'{index}.stderr'

    def success(out, err, operation):
        result = json.loads(out.read_bytes())
        require(err.stat().st_size == 0 and result['tool'] == PRODUCT and result['status'] == 'ok' and
                result['schema_version'] == '1.0' and result['tool_version'] == os.environ['PACKAGE_VERSION'] and
                result['result']['operation'] == operation, 'installed common success envelope differs')
        return result['result']

    out, _ = run(['--version'])
    require(out.read_text().strip() == PRODUCT + ' ' + os.environ['PACKAGE_VERSION'], 'installed version differs')
    out, _ = run(['--help'])
    require('normalize' in out.read_text() and '--summary' in out.read_text() and
            'summarize' not in out.read_text() and
            all(flag in out.read_text() for flag in ('--read-extension5', '--read-extension3', '--read2pos')),
            'installed counting/default-normalize scope differs')
    matrix = directory / 'counts.tsv'
    summary = directory / 'counts.summary.tsv'
    out, err = run(['library-a.sam', 'library-b.sam', '-a', 'genes.gtf', '-o', matrix,
                    '--summary', summary, '--threads', '1', '--json'])
    report = success(out, err, 'count')
    require(report['features'] == 3 and [p['total'] for p in report['libraries']] == [10, 9] and
            [p['assigned'] for p in report['libraries']] == [5, 4], 'installed count report differs')
    header = 'Geneid\tChr\tStart\tEnd\tStrand\tLength\tlibrary-a.sam\tlibrary-b.sam\n'
    body = ('geneA\tchr1;chr1\t101;201\t150;250\t+;+\t100\t3\t2\n'
            'geneB\tchr1\t141\t220\t-\t80\t1\t1\n'
            'geneC\tchr1;chr1\t301;341\t350;400\t+;+\t100\t1\t1\n')
    require(matrix.read_text() == '# Program:rsomics-count v' + os.environ['PACKAGE_VERSION'] + '\n' + header + body,
            'installed complete multi-library count matrix differs')
    categories = [('Assigned', 5), ('Unassigned_Unmapped', 1), ('Unassigned_Read_Type', 0),
                  ('Unassigned_Singleton', 0), ('Unassigned_MappingQuality', 0), ('Unassigned_Chimera', 0),
                  ('Unassigned_FragmentLength', 0), ('Unassigned_Duplicate', 0), ('Unassigned_MultiMapping', 1),
                  ('Unassigned_Secondary', 0), ('Unassigned_NonSplit', 0), ('Unassigned_NoFeatures', 1),
                  ('Unassigned_Overlapping_Length', 0), ('Unassigned_Ambiguity', 2)]
    require(summary.read_text() == 'Status\tlibrary-a.sam\tlibrary-b.sam\n' +
            ''.join(f'{name}\t{value}\t{value - (name == "Assigned")}\n' for name, value in categories),
            'installed complete assignment summary differs')
    projected_matrix = directory / 'projected.tsv'
    projected_summary = directory / 'projected.summary.tsv'
    out, err = run(['projection-paired.sam', '-a', 'projection-paired.saf', '-F', 'saf', '-o', projected_matrix,
                    '--summary', projected_summary, '--paired', '--strand', 'sense', '--multi-overlap', 'full',
                    '--read-extension5', '10', '--read-extension3', '7', '--read2pos', '5',
                    '--frac-overlap', '.05', '--threads', '1', '--json'])
    report = success(out, err, 'count')
    require(report['features'] == 2 and len(report['libraries']) == 1 and
            report['libraries'][0]['total'] == 1 and report['libraries'][0]['assigned'] == 1,
            'installed projected paired count report differs')
    require(projected_matrix.read_text() == '# Program:rsomics-count v' + os.environ['PACKAGE_VERSION'] + '\n' +
            'Geneid\tChr\tStart\tEnd\tStrand\tLength\tprojection-paired.sam\n'
            'five-forward\tchr1\t90\t100\t+\t11\t1\n'
            'five-reverse\tchr1\t130\t140\t+\t11\t1\n' and
            projected_summary.read_text() == 'Status\tprojection-paired.sam\n' +
            ''.join(f'{name}\t{int(name == "Assigned")}\n' for name, _ in categories),
            'installed complete projected paired matrix and summary differ')
    expected = {'tpm': [('571428.571429', '470588.235294'), ('238095.238095', '294117.647059'), ('190476.190476', '235294.117647')],
                'fpkm': [('6000000.000000', '5000000.000000'), ('2500000.000000', '3125000.000000'), ('2000000.000000', '2500000.000000')],
                'fpkm-uq': [('15000000.000000', '13333333.333333'), ('6250000.000000', '8333333.333333'), ('5000000.000000', '6666666.666667')]}
    for mode, values in expected.items():
        normalized = directory / (mode + '.tsv')
        out, err = run(['normalize', matrix, '-o', normalized, '--mode', mode, '--json'])
        report = success(out, err, 'normalize')
        require(report['genes'] == 3 and report['samples'] == 2 and report['mode'] == mode,
                'installed normalize report differs')
        require(normalized.read_text() == 'Geneid\tlibrary-a.sam\tlibrary-b.sam\n' +
                ''.join(f'{gene}\t{left}\t{right}\n' for gene, (left, right) in zip(('geneA', 'geneB', 'geneC'), values)),
                'installed complete normalized output differs: ' + mode)
    stdin_matrix = directory / 'stdin.tsv'
    out, err = run(['-', '-a', 'genes.gtf', '-o', stdin_matrix, '--threads', '1', '--json'],
                   stdin=fixture / 'single.sam')
    success(out, err, 'count')
    require(stdin_matrix.read_text() == '# Program:rsomics-count v' + os.environ['PACKAGE_VERSION'] + '\n' +
            header.replace('library-a.sam\tlibrary-b.sam', '-') +
            body.replace('\t3\t2\n', '\t3\n').replace('\t1\t1\n', '\t1\n') and
            Path(str(stdin_matrix) + '.summary').read_text() == 'Status\t-\n' +
            ''.join(f'{name}\t{value}\n' for name, value in categories), 'installed stdin/default summary differs')
    protected_matrix = directory / 'protected.tsv'
    protected_summary = directory / 'protected.summary.tsv'
    protected_matrix.write_bytes(b'previous counts\n')
    protected_summary.write_bytes(b'previous summary\n')
    malformed = directory / 'malformed.sam'
    malformed.write_bytes((fixture / 'single.sam').read_bytes() + b'broken\tSAM\n')
    out, err = run([malformed, '-a', 'genes.gtf', '-o', protected_matrix,
                    '--summary', protected_summary, '--json'], status=1)
    diagnostic = json.loads(err.read_bytes())
    require(out.stat().st_size == 0 and diagnostic['status'] == 'error' and diagnostic['tool'] == PRODUCT and
            diagnostic['error']['kind'] == 'InvalidInput' and diagnostic['exit_code'] == 1 and
            protected_matrix.read_bytes() == b'previous counts\n' and
            protected_summary.read_bytes() == b'previous summary\n',
            'installed malformed alignment lost diagnostics or paired output atomicity')
    malformed_counts = directory / 'malformed-counts.tsv'
    malformed_counts.write_bytes(b'Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample\ngene\tchr1\t1\t100\t+\t100\tbad\n')
    out, err = run(['normalize', malformed_counts, '-o', protected_matrix, '--json'], status=1)
    diagnostic = json.loads(err.read_bytes())
    require(out.stat().st_size == 0 and diagnostic['error']['kind'] == 'InvalidInput' and
            diagnostic['exit_code'] == 1 and protected_matrix.read_bytes() == b'previous counts\n',
            'installed malformed normalization lost diagnostics or prior output')
    require(len(records) == 10, 'installed lifecycle call count differs')
    save(directory / 'VERIFIED.json', dict(installed_binary=identity(binary), calls=records,
         default_count=True, projected_count=True, call_count=len(records),
         summary_side_output=True, normalize_modes=list(expected),
         own_output_reloaded=True, complete_matrix_and_summary=True, stdin_default_summary=True,
         common_json=True, malformed_alignment_preserves_both_outputs=True,
         malformed_normalization_preserves_output=True,
         files={p.name: identity(p) for p in directory.iterdir() if p.is_file()},
         scope='installed product lifecycle only; no local upstream oracle or performance execution'))

def bundle(binary, candidate=False):
    temporary = Path(os.environ['RUNNER_TEMP'])
    delivery = temporary / 'delivery'
    binary = Path(binary).resolve(strict=True)
    native_target = os.environ['NATIVE_TARGET']
    rustc = subprocess.check_output(['rustc', '-vV'], text=True)
    require(HOSTS[native_target] == os.environ['EXPECTED_HOST'] and
            '\nhost: ' + os.environ['EXPECTED_HOST'] + '\n' in rustc, 'bundle is not from the declared native host')
    stem = PRODUCT + '-' + os.environ['PACKAGE_VERSION'] + '-' + native_target
    folder = delivery / stem
    folder.mkdir()
    shutil.copy2(binary, folder / PRODUCT)
    for name in ('README.md', 'ATTRIBUTION.md', 'LICENSE-MIT', 'LICENSE-APACHE-2.0', 'THIRD_PARTY_NOTICES.md'):
        shutil.copy2(ROOT / name, folder / name)
    shutil.copytree(ROOT / 'licenses', folder / 'licenses')
    provenance = dict(product=PRODUCT, version=os.environ['PACKAGE_VERSION'], target=native_target,
         head_sha=os.environ['GITHUB_SHA'], run_id=int(os.environ['GITHUB_RUN_ID']),
         run_attempt=int(os.environ['GITHUB_RUN_ATTEMPT']), native_host=os.environ['EXPECTED_HOST'], rustc=rustc,
         machine=subprocess.check_output(['uname', '-a'], text=True).strip(), binary=identity(binary),
         lock_sha256=identity(ROOT / 'Cargo.lock')['sha256'],
         source_package=json.loads((delivery / 'package-verification.json').read_bytes()),
         native_dependencies=json.loads((delivery / 'native-dependencies.json').read_bytes()),
         installed_lifecycle=json.loads((delivery / 'lifecycle/VERIFIED.json').read_bytes()),
         accepted_implementation_base='2b5e450372689063736721dc1b807dcba9d62717',
         historical_performance_source='3a0154b15c8a985b999a98c4dd61d07a7533dec3',
         projection_performance_source='2b5e450372689063736721dc1b807dcba9d62717',
         projection_performance_record=dict(path='benchmarks/2026-11-03/projection-samples.json',
             **identity(ROOT / 'benchmarks/2026-11-03/projection-samples.json')),
         this_binary_remeasured=False, registry_publication=False, signed_release=False,
         immutable_registry_archive_equivalent=False, local_uncommitted_candidate=candidate)
    require(provenance['source_package']['clean_vcs'] == (not candidate) and provenance['source_package']['source_commit'] ==
            os.environ['GITHUB_SHA'], 'native bundle source package is not the clean CI head')
    require(provenance['installed_lifecycle']['installed_binary'] == identity(binary), 'bundle/lifecycle binary differs')
    save(folder / 'BUILD.json', provenance)
    paths = [folder, *folder.rglob('*')]
    require(len(paths) == len({str(p.relative_to(folder)).casefold() for p in paths}), 'case-colliding bundle path')
    for path in paths:
        require(not path.is_symlink() and (path.is_dir() or path.is_file()), 'unsupported bundle member')
        path.chmod(0o755 if path.is_dir() or path.name == PRODUCT else 0o644)
    archive = delivery / (stem + '.tar.gz')
    with tarfile.open(archive, 'w:gz') as output:
        output.add(folder, arcname=stem)
    source_archive = Path(os.environ['SOURCE_PACKAGE'])
    (delivery / 'SHA256SUMS').write_text(''.join(identity(path)['sha256'] + '  ' + path.name + '\n'
                                              for path in (archive, source_archive)))
    print(json.dumps(provenance, indent=2, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('source', 'smoke', 'bundle'))
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--candidate', action='store_true', help='explicitly uncommitted local packaging; never a release gate')
    args = parser.parse_args()
    if args.action == 'source':
        source(args.candidate)
    elif args.action == 'smoke':
        smoke(args.binary)
    else:
        bundle(args.binary, args.candidate)


if __name__ == '__main__':
    main()
