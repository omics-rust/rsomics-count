#!/bin/bash
set -euo pipefail
/usr/bin/python3 -B - <<'PY'
import os
import re
import sys

overrides = {
    'CARGO_BUILD_TARGET', 'CARGO_BUILD_BUILD_DIR', 'CARGO_BUILD_RUSTC_WRAPPER',
    'CARGO_BUILD_RUSTC_WORKSPACE_WRAPPER', 'CARGO_BUILD_RUSTFLAGS',
    'RUSTC_WRAPPER', 'RUSTC_WORKSPACE_WRAPPER', 'RUSTFLAGS', 'CARGO_ENCODED_RUSTFLAGS',
    'RUSTDOCFLAGS', 'CARGO_ENCODED_RUSTDOCFLAGS', 'CC', 'CXX', 'CFLAGS', 'CXXFLAGS',
    'CPPFLAGS', 'CXXSTDLIB', 'AR', 'ARFLAGS', 'RANLIB', 'RANLIBFLAGS', 'CRATE_CC_NO_DEFAULTS',
}
for name, value in os.environ.items():
    if value and (name in overrides or name.startswith((
        'CC_', 'CXX_', 'CFLAGS_', 'CXXFLAGS_', 'CPPFLAGS_', 'CXXSTDLIB_',
        'AR_', 'ARFLAGS_', 'RANLIB_', 'RANLIBFLAGS_')) or re.fullmatch(
        r'(HOST|TARGET)_(CC|CXX|CFLAGS|CXXFLAGS|CPPFLAGS|CXXSTDLIB|AR|ARFLAGS|RANLIB|RANLIBFLAGS)', name)
        or re.fullmatch(r'CARGO_TARGET_.+_(RUSTFLAGS|LINKER|RUNNER)', name)):
        sys.exit('Unsupported compiler/target override: ' + name)
PY
export CARGO_HOME=/Volumes/KIOXIA/Developments/cargo-home
export CARGO_TARGET_DIR=/Volumes/KIOXIA/Developments/cargo-target/rsomics-count
export RUSTUP_HOME=/Volumes/KIOXIA/Developments/rustup-home
export TMPDIR=/Volumes/KIOXIA/Developments/tmp
export XDG_CACHE_HOME=/Volumes/KIOXIA/Developments/tmp/cache
export GH_CACHE_DIR=/Volumes/KIOXIA/Developments/tmp/gh-cache
export CLANG_MODULE_CACHE_PATH=/Volumes/KIOXIA/Developments/tmp/clang-modules-count
export CCACHE_DIR=/Volumes/KIOXIA/Developments/tmp/ccache-count
export SCCACHE_DIR=/Volumes/KIOXIA/Developments/tmp/sccache-count
export RUSTC_WRAPPER= RUSTC_WORKSPACE_WRAPPER=
export CARGO_ENCODED_RUSTFLAGS= CARGO_ENCODED_RUSTDOCFLAGS=
export PYTHONDONTWRITEBYTECODE=1
export CARGO_BUILD_JOBS=2
export CARGO_INCREMENTAL=0
toolchain="$RUSTUP_HOME/toolchains/1.91.0-aarch64-apple-darwin/bin"
export CARGO="$toolchain/cargo"
export RUSTC="$toolchain/rustc"
export RUSTDOC="$toolchain/rustdoc"
export PATH="$toolchain:$PATH"
for argument in "$@"; do
  case "$argument" in --target-dir|--target-dir=*|--config|--config=*) echo "Storage overrides are not supported: $argument" >&2; exit 1;; esac
done
/usr/sbin/diskutil info -plist /Volumes/KIOXIA | /usr/bin/python3 -B -c '
import plistlib, sys
volume = plistlib.loads(sys.stdin.buffer.read())
if volume.get("MountPoint") != "/Volumes/KIOXIA" or volume.get("Internal") is not False:
    sys.exit("KIOXIA is not mounted as an external volume")
print("verified external volume:", volume["DeviceNode"], file=sys.stderr)
'
external_device=$(stat -f %d /Volumes/KIOXIA)
for path in "$CARGO_HOME" "$RUSTUP_HOME" "$TMPDIR" "$toolchain" /Volumes/KIOXIA/Developments/cargo-target; do
  resolved=$(cd "$path" && pwd -P)
  case "$resolved" in /Volumes/KIOXIA/Developments/*) ;; *) echo "Forbidden build path: $resolved" >&2; exit 1;; esac
  test "$(stat -f %d "$resolved")" = "$external_device" || { echo "Path is not on the verified volume: $resolved" >&2; exit 1; }
done
available_kib=$(df -Pk /Volumes/KIOXIA | awk 'NR == 2 {print $4}')
test "$available_kib" -ge 8388608 || { echo "External disk has less than 8 GiB independent headroom" >&2; exit 1; }
for path in "$CARGO_TARGET_DIR" "$XDG_CACHE_HOME" "$GH_CACHE_DIR" \
    "$CLANG_MODULE_CACHE_PATH" "$CCACHE_DIR" "$SCCACHE_DIR"; do
  if [[ ! -e "$path" && ! -L "$path" ]]; then mkdir "$path"; fi
  resolved=$(cd "$path" && pwd -P)
  case "$resolved" in /Volumes/KIOXIA/Developments/*) ;; *) echo "Forbidden output/cache path: $resolved" >&2; exit 1;; esac
  test "$(stat -f %d "$resolved")" = "$external_device" || { echo "Output/cache is not on the verified volume: $resolved" >&2; exit 1; }
  echo "verified build/cache path: $resolved" >&2
done
for path in "$CARGO" "$RUSTC" "$RUSTDOC" "$toolchain/cargo-clippy"; do
  test -x "$path" && test "$(stat -f %d "$path")" = "$external_device" || { echo "Tool is not on the verified volume: $path" >&2; exit 1; }
done
df -h /Volumes/KIOXIA >&2
target=$("$CARGO" metadata --locked --offline --no-deps --format-version 1 | /usr/bin/python3 -B -c 'import json,sys; print(json.load(sys.stdin)["target_directory"])')
test "$target" = "$CARGO_TARGET_DIR"
echo "Cargo resolved target: $target" >&2
if [[ "${1:-}" == "--exec" ]]; then shift; exec "$@"; fi
if [[ "${1:-}" == "clippy" ]]; then shift; exec "$toolchain/cargo-clippy" clippy "$@"; fi
exec "$CARGO" "$@"
