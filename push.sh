#!/bin/bash

TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

git add .
git commit -m "$TIMESTAMP"
git branch -M main
git remote add origin git@github.com:larena14/NS-Attack1Scenario.git
git push -u origin main
