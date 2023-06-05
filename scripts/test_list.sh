#!/bin/bash

cc='contour_'
outmap='_labelmap'
border='_border'

while read line  
do   
   I="${line%.*}"
   ./SCALP -i $2$line -k $3 -m $4 -outm $I$outmap.png -outb $I$border.png -c $2$cc$I.png \n
done < $1
