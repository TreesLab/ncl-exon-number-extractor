### Examples

```
cat NCLscan.result | python ../ncl_exon_number_extractor.py gencode.v28.annotation.gtf > out.default.tsv
cat NCLscan.result | python ../ncl_exon_number_extractor.py gencode.v28.annotation.gtf --show-all > out.show_all.tsv
cat NCLscan.result | python ../ncl_exon_number_extractor.py gencode.v28.annotation.gtf --show-all --expand > out.show_all_expand.tsv

cat NCLscan.result | python ../ncl_exon_number_extractor.py gencode.v28.annotation.gtf --detail > out.detail.tsv
cat NCLscan.result | python ../ncl_exon_number_extractor.py gencode.v28.annotation.gtf --show-all --detail > out.show_all_detail.tsv
cat NCLscan.result | python ../ncl_exon_number_extractor.py gencode.v28.annotation.gtf --show-all --expand --detail > out.show_all_expand_detail.tsv
```

