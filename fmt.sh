#!/bin/bash

FILES=()
while IFS=  read -r -d $'\0'; do
    FILES+=("$REPLY")
done < <(sh -c "find . -regex '.*\.\(cpp\|hpp\|cc\|cxx\|h\)' -print0")

for F in ${FILES[@]}; do
    echo "$F"
    clang-format -i $F
done

echo "done"
