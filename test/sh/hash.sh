#!/bin/bash
declare -A animals
animals=( ["moo"]="cow" ["bird"]="chirp" )

echo "${animals["moo"]}"
for sound in "${!animals[@]}"; do echo "$sound - ${animals["$sound"]}"; done

