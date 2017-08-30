# Pipeline para desenvolvimento de sequências para Phage Display

Para o desenvolvimento de um pipeline de sequências das sequencias a serem utilizadas no experimento de phage display, eu criei este repositório que contém o código de preenchimento e trimming de sequências.

O pipeline inteiro é:
```bash
pyfasta split -n 1 -k 156 -o 21 FASTAFILE.fasta
python phagedisplay/fasta_normalizer.py -f FASTAFILE_OUTPUT_FROM_PYFASTA_HERE.fasta -l 156 -o 21 -of FASTAFILE-OUTPUTNAME_HERE.fasta

```

## fasta_normalizer.py

Para utilizar o fasta_normalizer.py na linha de comando, você deve passar os seguintes parâmetros:

```
Você deve escolher uma ação que deve ser falada no inicio do comando
-n ou --fasta_normalizer - Para normalizar o fasta e
-p ou --peptide-fasta - Transformar um fasta de DNA para Peptídeos
-c ou --fasta-cleaner - Para realizar a limpeza do fasta (remoção das sequências com alfabeto estranho)

E preencher com os dados que são requeridos por ela
-f ou -file : Arquivo de input
-l ou -length : Tamanho desejado das sequencias
-o ou -overlap : Tamanho do overlap das sequencias
-of ou --output_file : Arquivo de output das sequencias
```
