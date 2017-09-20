# Pipeline para desenvolvimento de sequências para Phage Display

Para o desenvolvimento de um pipeline de sequências das sequencias a serem utilizadas no experimento de phage display, eu criei este repositório que contém o código de preenchimento e trimming de sequências.

O pipeline inteiro é:
```bash
sed "s/\*$//" ORIGINAL-FASTAFILE.FASTA > PRE_FASTAFILE.fasta
sed "s/|.*//" ORIGINAL-FASTAFILE.FASTA > PRE_FASTAFILE.fasta
pyfasta split -n 1 -k 52 -o 7 FASTAFILE.fasta
sed "s/\([0-9]\)_\(-\?[0-9]\)/\1\|\2/g" PRE_FASTAFILE.fasta > FASTAFILE.fasta
python phagedisplay/fasta_normalizer.py -n -f FASTAFILE_OUTPUT_FROM_PYFASTA_HERE.fasta -l 52 -o 7 -of FASTAFILE-OUTPUTNAME_HERE.fasta
python phagedisplay/fasta_normalizer.py -c -f FASTAFILE-OUTPUTNAME_HERE.fasta -w x -of FASTAFILE-OUTPUTNAME_HERE.fasta.clean
python phagedisplay/fasta_normalizer.py -p -f FASTAFILE-OUTPUTNAME_HERE.fasta.clean -of FASTAFILE-OUTPUTNAME_HERE.fasta.pep
makeblastdb -in FASTAFILE-OUTPUTNAME_HERE.fasta.pep -dbtype prot
blastp -task blastp -db FASTAFILE-OUTPUTNAME_HERE.fasta.pep -query FASTAFILE-OUTPUTNAME_HERE.fasta.pep -outfmt 6 -qcov_hsp_perc 100 -out blast-peptide.blast
awk '{if ($3 == 100) print $0}' blast-peptide.blast > blastOUT-filtered100id.blast
```

## fasta_normalizer.py

Para utilizar o fasta_normalizer.py na linha de comando, você deve passar os seguintes parâmetros:

```
Você deve escolher uma ação que deve ser falada no inicio do comando
-n ou --fasta_normalizer - Para normalizar o fasta para o mesmo tamanho de sequencias
-p ou --peptide-fasta - Transformar um fasta de DNA para Peptídeos
-c ou --fasta-cleaner - Para realizar a limpeza do fasta (remoção das sequências com alfabeto estranho)

E preencher com os dados que são requeridos por ela
-f ou -file : Arquivo de input
-l ou -length : Tamanho desejado das sequencias
-o ou -overlap : Tamanho do overlap das sequencias
-of ou --output_file : Arquivo de output das sequencias
```
