#!/bin/bash

file='DB/zinc_1k_50k/zinc100k.smi' #현재 .sh 파일 위치가 /workspace/joonseon 일 경우에만 해당함.

save_path='DB/zinc_1k_50k/pdbqt'

#mkdir -p pdbqt/$dir #-p 옵션으로 pdbqt 하위에 여러 하위 디렉토리 지정해서 만듦
files=$(ls $dir)  # file load, $(....) 변수로 (...)의 명령이 실행되게 변경

cnt=0

while IFS= read -r line
do
  obabel -:$line -o mol2 -h --gen3D -O $save_path/mol_$cnt.mol2 #$dir/을 추가하여 smi_file 의 위치를 정확하게 바꿔줌. $smi_file 로만 해 주면 디렉토리가 없어 파일 없으므로 실행 안 됨.
  obabel -i mol2 $save_path/mol_$cnt.mol2 -o pdbqt -O $save_path/mol_$cnt.pdbqt
  echo $cnt
  ((cnt++))
done < $file


dir='DB/zinc_1k_50k/' #현재 .sh 파일 위치가 /workspace/joonseon 일 경우에만 해당함.
# save_path='$dir/pdbqt'

files=$(ls $dir/pdbqt)

for pdbqt_file in ${files[@]};do
   ID=$(echo $pdbqt_file | cut -d '/' -f 2 | cut -d '.' -f 1) #DB/ID.smi -> ID.smi 로 출력 => ID.smi를 ID로 자르기; 두 줄 명령을 한 줄로 짧게 변경
   ./vina --config ./target/5rlz_dock.cfg --ligand ./DB/zinc_1k_50k/pdbqt/$ID.pdbqt --out ./DB/zinc_1k_50k/docking_result/${ID}_log.txt

   done
