#!/usr/bin/env bash

if [ "$SAGE_LOCAL" = "" ]; then
   echo "SAGE_LOCAL undefined ... exiting";
   echo "Maybe run 'sage -sh'?"
   exit 1
fi

cd src

./configure --prefix="$SAGE_LOCAL"
if [ $? -ne 0 ]; then
   echo "Error configuring PACKAGE_NAME."
   exit 1
fi

python3 -m build