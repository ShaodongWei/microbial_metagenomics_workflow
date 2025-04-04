mamba has to be installed also, due to newer version of snakemake 
use '--use-conda --conda-frontend conda', if mamba is problematic to run, e.g. parameter not recognized, then choose to use conda. so far I cannot make mamba work. 

use channel_priority: strict inside yaml

using mamba will lead to --no-default-packages not recognized error. 
This is because Mamba 2.0.0, released last week, no longer supports the --no-default-packages option that Snakemake passes when creating an environment.

we should use snakemake binning --cores 40 --use-conda --conda-frontend conda --resources gpu=2, so that binning can have GPU
