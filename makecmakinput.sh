#!/bin/bash
cat $1 | grep -v "#" | cut -d "=" -f 2 | sed -e  "s/\\\\/ /g" -e "s/\.o/\.cpp;/g" | tr "\n" " " | tr -d "\t\n\r " | tr ";" " " > $2