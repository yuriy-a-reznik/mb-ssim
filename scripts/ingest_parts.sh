#!/bin/sh
# Assemble a database that arrived as multiple uploaded archive parts.
#
#   scripts/ingest_parts.sh SRCDIR DESTDIR
#
# Handles the three shapes these usually arrive in:
#   * split zip volumes  (name.z01 name.z02 ... name.zip)  -> recombined first
#   * independent zips   (refs.zip dist1.zip dist2.zip)    -> unpacked in place
#   * a single rar       (tid2013.rar)                     -> unrar / 7z
#
# Everything is extracted into DESTDIR, then the layout is reported so the
# right make_index.py invocation is obvious.

set -e
SRC="${1:?usage: $0 SRCDIR DESTDIR}"
DST="${2:?usage: $0 SRCDIR DESTDIR}"
mkdir -p "$DST"

echo "== source archives"
ls -la "$SRC" | grep -Ei '\.(zip|z[0-9]+|rar|r[0-9]+|7z|tar|gz)$' || echo "(none found)"

# Split zip volumes: recombine before unpacking.
for main in "$SRC"/*.zip; do
    [ -e "$main" ] || continue
    base="${main%.zip}"
    if [ -e "${base}.z01" ]; then
        echo "== recombining split volumes for $(basename "$base")"
        zip -q -s 0 "$main" --out "$DST/_combined.zip"
        unzip -q -o "$DST/_combined.zip" -d "$DST"
        rm -f "$DST/_combined.zip"
        touch "$DST/.handled_$(basename "$base")"
    fi
done

# Plain zips (skip any that were part of a split set).
for z in "$SRC"/*.zip; do
    [ -e "$z" ] || continue
    base="${z%.zip}"
    [ -e "$DST/.handled_$(basename "$base")" ] && continue
    echo "== unzip $(basename "$z")"
    unzip -q -o "$z" -d "$DST"
done

# RAR.
for r in "$SRC"/*.rar; do
    [ -e "$r" ] || continue
    echo "== unrar $(basename "$r")"
    if command -v unrar >/dev/null; then unrar x -inul -o+ "$r" "$DST/"
    else 7z x -bso0 -y -o"$DST" "$r"; fi
done

rm -f "$DST"/.handled_*

echo
echo "== extracted layout (two levels)"
find "$DST" -maxdepth 2 -type d | head -40
echo
echo "== file counts by extension"
find "$DST" -type f | sed 's/.*\.//' | tr 'A-Z' 'a-z' | sort | uniq -c | sort -rn | head
echo
echo "== score files present"
find "$DST" -maxdepth 3 -type f \( -iname 'mos*' -o -iname '*.mat' -o -iname '*.xlsx' \) | head
