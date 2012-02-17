A="1 2 3   4     7  "
B=$(echo $A | sed -e 's/[ ]/,/g')
echo $B 
