#!/bin/bash
rm exam/exe
rm exam/*.out
rm exam/*.o
cd exam
make
rm *.o
cd ..
rm -rf data/*
rm -rf Tem*
for a in {1..21} 
do
mkdir Tem$a
cp -r exam/exe Tem$a
mkdir Tem$a/buffer
echo $a >Tem$a/m1.dat
echo -e "#!/bin/bash\ncd Tem$a\n./exe" >Tem$a/run.sh
chmod 755 run.sh
done
sh Tem1/run.sh &
sh Tem2/run.sh &
sh Tem3/run.sh &
sh Tem4/run.sh &
sh Tem5/run.sh &
sh Tem6/run.sh &
sh Tem7/run.sh &
sh Tem8/run.sh &
sh Tem9/run.sh &
sh Tem10/run.sh &
sh Tem11/run.sh &
sh Tem12/run.sh &
sh Tem13/run.sh &
sh Tem14/run.sh &
sh Tem15/run.sh &
sh Tem16/run.sh &
sh Tem17/run.sh &
sh Tem18/run.sh &
sh Tem19/run.sh &
sh Tem20/run.sh &
sh Tem21/run.sh &
wait
for a in {1..21}
do
mv Tem$a/buffer/Vtotal.dat data/V$a.dat
done
cd final
rm *.o
rm *.out
rm exe
cd buffer
rm *.dat
cd ..
make
rm *.o
./exe