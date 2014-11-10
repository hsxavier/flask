#!/bin/bash
# This script runs gencovl e em seguida corrlnfields. 
# Ambos dependem do mesmo arquivo de configuração que é passado como argumento.

echo "../bin/gencovl $1"
../bin/gencovl $1

echo "../bin/corrlnfields $1"
../bin/corrlnfields $1

