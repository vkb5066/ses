#!/usr/bin/csh

set k=2

while($k < 26)
	echo "$k"

	echo "com" >kpoints
	echo "m" >>kpoints	
	echo "$k $k $k" >>kpoints

	./tst >/dev/null
	./calcbg.py

	@ k++
	@ k++
end
