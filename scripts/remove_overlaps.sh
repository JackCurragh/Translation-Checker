

for file in cls_annotations/translation/*.overlapping.bed; do
    BASE=$( basename $file .overlapping.bed )
    echo $BASE
    cp $file ./cls_annotations_bak/${BASE}.summary.csv
    cat $file | while read line; do 

        grep -v "$line" cls_annotations/translation/$BASE.summary.csv > cls_annotations/translation/$BASE.summary.csv.tmp
        mv cls_annotations/translation/$BASE.summary.csv.tmp cls_annotations/translation/$BASE.summary.csv
    done
done