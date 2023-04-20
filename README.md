# LigandDiversity

Small introduction to generating a diversified subset of ligands and decoys (in this example: taken from the DUD-E)
If you want to use the script for both, please change the files names so they aren't overwritten.
Step by Step:

1. Download Anaconda / Miniconda (Miniconda in this example)
2. Open the Miniconda command prompt from your search bar or .exe (Anaconda Prompt (Miniconda 3))
3. Type in:
	conda create -n **NAME** -c conda-forge python jupyter jupyterlab rdkit

	-> hit enter
	-> this creates a new working environment and installs jupyter, jupyterlab, python and rdkit into that specific environment
4. You are now in the environment (when closing and opening again THROUGH THE ANACONDA PROMPT, type in: "activate **NAME**" to enter your environment again
5. When in your environment, typing in 'jupyter' opens up jupyter-notebook, using 'jupyter-lab' opens jupyterlab -> I recommend jupyterlab, as it helps with overview
6. copy & paste the code and it should be able to generate your files -> change your paths accordingly, change names for intermediate files as you wish -> not every file has to be written, here it's for completion's sake
7. Watch the indentations!
8. In the last step - diversity - you can change your needed sample size accordingly; default is 30 compounds "line: lazypick" change the number to the desired number
9. You can run each cell seperately with Shift + Enter; writing all the blocks in different cells helps trouble shooting as you can run each step separately. 
10. Have fun :)
