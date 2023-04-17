

function backup_and_remove_lines() {
    local file="$1"
    local BASE="$( basename "$file" .overlapping.bed)"

    cp "$file" "./cls_annotations_bak/${BASE}.summary.csv"

    while read -r line; do
        grep -v "$line" "cls_annotations/translation/$BASE.summary.csv" > "cls_annotations/translation/$BASE.summary.csv.tmp"
        mv "cls_annotations/translation/$BASE.summary.csv.tmp" "cls_annotations/translation/$BASE.summary.csv"
    done < "$file"
}

for file in cls_annotations/translation/*.overlapping.bed; do
    BASE=$(basename "$file" .overlapping.bed)
    # echo "$(wc -l "$file") overlaps found in $BASE)"
    backup_and_remove_lines "$file" &
    # echo "$(wc -l "cls_annotations/translation/$BASE.summary.csv") lines remaining in $BASE)"
done