#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
REPO_ROOT="$(CDPATH= cd -- "$SCRIPT_DIR/../.." && pwd)"

BINARY_PATH="${1:-$REPO_ROOT/src/COMPAS}"
ARCH_NAME="${2:-$(uname -m)}"
BUNDLE_DIR="${3:-$REPO_ROOT/dist/COMPAS-macos-$ARCH_NAME}"
BIN_DIR="$BUNDLE_DIR/bin"
LIB_DIR="$BUNDLE_DIR/lib"
README_PATH="$BUNDLE_DIR/README.txt"
SEEN_DEPS="$(mktemp)"

cleanup() {
    rm -f "$SEEN_DEPS"
}

trap cleanup EXIT

if [ ! -x "$BINARY_PATH" ]; then
    echo "Expected executable COMPAS binary at '$BINARY_PATH'" >&2
    exit 1
fi

for required_command in codesign install_name_tool otool; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "$required_command is required to package the macOS bundle" >&2
        exit 1
    fi
done

maybe_strip() {
    local path="$1"

    if ! command -v strip >/dev/null 2>&1; then
        return 0
    fi

    strip -x "$path" 2>/dev/null || true
}

is_system_library() {
    case "$1" in
        /System/*|/usr/lib/*)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

dependencies_for() {
    otool -L "$1" | tail -n +2 | awk '{print $1}'
}

copy_dependency() {
    local source_path="$1"
    local dest_path="$LIB_DIR/$(basename "$source_path")"

    if [ ! -f "$dest_path" ]; then
        cp -f "$source_path" "$dest_path"
        chmod 755 "$dest_path"
        maybe_strip "$dest_path"
    fi
}

rewrite_references() {
    local target="$1"
    local target_dir_mode="$2"
    local dependency
    local dependency_name

    while IFS= read -r dependency; do
        dependency_name="$(basename "$dependency")"
        if [ -f "$LIB_DIR/$dependency_name" ]; then
            if [ "$target_dir_mode" = "binary" ]; then
                install_name_tool -change "$dependency" "@loader_path/../lib/$dependency_name" "$target" || true
            else
                install_name_tool -change "$dependency" "@loader_path/$dependency_name" "$target" || true
            fi
        fi
    done < <(dependencies_for "$target")
}

rm -rf "$BUNDLE_DIR"
mkdir -p "$BIN_DIR" "$LIB_DIR"

cp -f "$BINARY_PATH" "$BIN_DIR/COMPAS"
cp "$SCRIPT_DIR/run_compas.sh" "$BUNDLE_DIR/run_compas.sh"
chmod 755 "$BIN_DIR/COMPAS" "$BUNDLE_DIR/run_compas.sh"
maybe_strip "$BIN_DIR/COMPAS"

declare -a queue=("$BIN_DIR/COMPAS")

while [ "${#queue[@]}" -gt 0 ]; do
    current_file="${queue[0]}"
    queue=("${queue[@]:1}")

    while IFS= read -r dependency; do
        case "$dependency" in
            @loader_path/*|@executable_path/*|@rpath/*|"")
                continue
                ;;
        esac

        if is_system_library "$dependency"; then
            continue
        fi

        dependency_name="$(basename "$dependency")"
        copy_dependency "$dependency"

        if ! grep -qx "$dependency_name" "$SEEN_DEPS" 2>/dev/null; then
            echo "$dependency_name" >> "$SEEN_DEPS"
            queue+=("$LIB_DIR/$dependency_name")
        fi
    done < <(dependencies_for "$current_file")
done

for library_path in "$LIB_DIR"/*; do
    [ -f "$library_path" ] || continue
    library_name="$(basename "$library_path")"
    install_name_tool -id "@loader_path/$library_name" "$library_path" || true
    rewrite_references "$library_path" "library"
    codesign --force --sign - "$library_path" >/dev/null 2>&1 || true
done

rewrite_references "$BIN_DIR/COMPAS" "binary"
codesign --force --sign - "$BIN_DIR/COMPAS" >/dev/null 2>&1 || true

cat > "$README_PATH" <<EOF
COMPAS bundled macOS artifact

Supported platform:
- macOS
- Architecture: $ARCH_NAME

Contents:
- bin/COMPAS: bundled COMPAS executable
- lib/: shared libraries packaged with this build
- run_compas.sh: launcher for the bundled executable

How to run:
1. Unpack this directory on a supported macOS machine.
2. From inside the unpacked directory, run:
   ./run_compas.sh -v

Notes and caveats:
- This bundle is intended to run without separately installing Boost, GSL, or HDF5.
- If Gatekeeper or local signing policy blocks execution, inspect the ad-hoc signed files in bin/ and lib/.
EOF

echo
echo "Created bundle at: $BUNDLE_DIR"
echo "Bundled shared libraries:"
find "$LIB_DIR" -maxdepth 1 -type f | sort
