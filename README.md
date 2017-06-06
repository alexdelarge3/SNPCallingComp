# SNPCallingComp

To make your graph:
snakemake --dag --configfile config.json -s workflow/SnakeBromance.mk | dot -Tpdf > dag.pdf

to run:
snakemake --configfile config.json -s workflow/SnakeBromance.mk
