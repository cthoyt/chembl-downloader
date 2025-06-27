# sort the bib file
biber --tool --output_align --output_indent=2 --output_fieldcase=lower paper.bib -O paper.bib

# see the instructions at
# https://joss.readthedocs.io/en/latest/submitting.html#docker
docker run --rm \
  --volume $PWD:/data \
  --user $(id -u):$(id -g) \
  --env JOURNAL=joss \
  openjournals/inara \
  -o pdf,preprint \
  paper.md

latexmk -pdf paper.preprint.tex
rm paper.bib.blg
rm paper.preprint.aux
rm paper.preprint.fdb_latexmk
rm paper.preprint.fls
rm paper.preprint.log
rm paper.preprint.tex
