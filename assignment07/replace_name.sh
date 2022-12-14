
# What to replace from
ORIG="YOUR_NAME"
# Names to replace with
NAME1="GIESEL"
NAME2="MUEHLBAUER"


# REPLACE FILENAME
for f in *.java; do
	mv "$f" "$(echo "$f" | sed s/$ORIG/$NAME1\_$NAME2/)"; 
done


# REPLACE INSIDE FILE
for f in *.java; do 
	sed -i "s/YOUR_NAME/$NAME1\_$NAME2/g" $f
done
