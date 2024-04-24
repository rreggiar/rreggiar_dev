# Analysis tools to process Parse Bioscience's Single Cell data.

Version v1.2.1c

## For detailed installation and use

See company page for background information, usage instructions and license.
https://www.parsebiosciences.com/

## Minimal install instructions:

Activate conda

Run install scripts for dependencies:

`bash ./install_dependencies_conda.sh -i -y`

IF you are doing TCR analysis:

`bash ./install_TCR.sh -i -y`

Install the pipeline per se:

`pip install . --no-cache-dir`

Test split-pipe (usage screen):

`split-pipe -h`

