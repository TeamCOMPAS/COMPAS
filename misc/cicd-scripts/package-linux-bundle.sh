#!/usr/bin/env bash

# Create the redistributable Linux COMPAS bundle used by the native tarball
# workflow and by the Linux PyPI wheel build.

set -euo pipefail

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
REPO_ROOT="$(CDPATH= cd -- "$SCRIPT_DIR/../.." && pwd)"

BINARY_PATH="${1:-$REPO_ROOT/src/COMPAS}"
BUNDLE_DIR="${2:-$REPO_ROOT/dist/COMPAS-linux-x86_64}"
BIN_DIR="$BUNDLE_DIR/bin"
LIB_DIR="$BUNDLE_DIR/lib"
README_PATH="$BUNDLE_DIR/README.txt"
TMP_LDD="$(mktemp)"

cleanup() {
    rm -f "$TMP_LDD"
}

trap cleanup EXIT

maybe_strip() {
    local path="$1"

    if ! command -v strip >/dev/null 2>&1; then
        return 0
    fi

    # Best-effort size reduction only. Some binaries or shared libraries may
    # already be stripped or may not support the requested mode.
    strip --strip-unneeded "$path" 2>/dev/null || strip "$path" 2>/dev/null || true
}

if [ ! -x "$BINARY_PATH" ]; then
    echo "Expected executable COMPAS binary at '$BINARY_PATH'" >&2
    exit 1
fi

if ! command -v ldd >/dev/null 2>&1; then
    echo "ldd is required to package the Linux bundle" >&2
    exit 1
fi

rm -rf "$BUNDLE_DIR"
mkdir -p "$BIN_DIR" "$LIB_DIR"

cp -Lf "$BINARY_PATH" "$BIN_DIR/COMPAS"
cp "$SCRIPT_DIR/run_compas.sh" "$BUNDLE_DIR/run_compas.sh"
chmod 755 "$BIN_DIR/COMPAS" "$BUNDLE_DIR/run_compas.sh"
maybe_strip "$BIN_DIR/COMPAS"

ldd "$BIN_DIR/COMPAS" | tee "$TMP_LDD"

if grep -q 'not found' "$TMP_LDD"; then
    echo "Cannot package bundle with unresolved runtime dependencies" >&2
    exit 1
fi

awk '
    /=>/ && $3 ~ /^\// { print $3; next }
    $1 ~ /^\// { print $1; next }
' "$TMP_LDD" | sort -u | while IFS= read -r lib; do
    case "$lib" in
        /lib64/ld-linux-*|/lib/x86_64-linux-gnu/ld-linux-*|*/libc.so.*|*/libm.so.*|*/libpthread.so.*|*/libdl.so.*|*/librt.so.*|*/libresolv.so.*|*/libnss_*.so.*|*/libcrypt.so.*|*/libutil.so.*|*/libanl.so.*)
            echo "Skipping system runtime: $lib"
            ;;
        *)
            cp -Ln "$lib" "$LIB_DIR/"
            maybe_strip "$LIB_DIR/$(basename "$lib")"
            ;;
    esac
done

cat > "$README_PATH" <<'EOF'
COMPAS bundled Linux artifact

Supported platform:
- Ubuntu 22.04
- Linux x86_64

Contents:
- bin/COMPAS: bundled COMPAS executable
- lib/: shared libraries packaged with this build
- run_compas.sh: launcher that sets LD_LIBRARY_PATH for the bundled libs

How to run:
1. Unpack this directory on a supported Linux x86_64 machine.
2. From inside the unpacked directory, run:
   ./run_compas.sh -v

Notes and caveats:
- This bundle is intended to run without separately installing Boost, GSL, or HDF5.
- It still relies on the target machine's glibc and other base system libraries being compatible with Ubuntu 22.04.
- If you invoke bin/COMPAS directly, the bundled lib/ directory will not be added automatically; use run_compas.sh instead.
EOF

echo
echo "Created bundle at: $BUNDLE_DIR"
echo "Bundled shared libraries:"
find "$LIB_DIR" -maxdepth 1 -type f | sort
