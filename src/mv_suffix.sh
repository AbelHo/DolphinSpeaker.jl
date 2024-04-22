#!/bin/bash

src="/Volumes/One Touch/res/concretecho/outvid"
dst="/Volumes/One Touch/res/concretecho/outvid_timeplot"
suffix="*_timeplot.mp4"
a="*_vid-combined.mp4"

find "$src" -type f -name $suffix -print0 | while IFS= read -r -d '' file; do
    dir=$(dirname "${file#$src}")
    mkdir -p "$dst/$dir"
    mv "$file" "$dst/$dir"
done