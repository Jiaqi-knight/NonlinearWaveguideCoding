#!/bin/bash
cd dados/
#for i in $(seq 1 1 15000)
#do
	#var="$i.mat"
	#rm $var
   	#echo "Welcome $i times"
#done

for i in $(seq 15001 1 80000)
do
	nome_antigo="$i.mat"
	nome_novo="$((i-15000)).mat"
	mv $nome_antigo $nome_novo
   	echo "Renomeando $nome_novo arquivo"
done