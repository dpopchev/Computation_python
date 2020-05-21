# set bash as default shell
SHELL := /usr/bin/env bash

# set default python
PYTHON 	 := /usr/bin/env python3
UNITTEST := $(PYTHON) -m unittest -v

# ln should always overwrite the created symbolic links
LN := ln
LN += --force
LN += --symbolic

MKDIR := mkdir
MKDIR += --parents

RM := rm
RM += --force
RM += --recursive

.PHONY: venv run clean

# make $(file) the default target for testing
# file is passed as argument, i.e. make default file=your_file_of_choice
venv:
	@$(MKDIR) bin/data
	@unlink bin/main.py
	@$(LN) ../src/$(file) bin/main.py

run:
	$(PYTHON) bin/main.py

clean:
	@$(RM) bin
