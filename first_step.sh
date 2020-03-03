#!/bin/bash

period="201204"


#rsync -a --progress /paula1/paula*/$period*   /paula3/QUEST/data
#rsync -a --progress /paula2/rac72/$period*   /paula3/QUEST/data

gunzip -r /paula3/QUEST/data/$period*

