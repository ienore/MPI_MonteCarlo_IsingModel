## MPI_MonteCarlo_IsingModel
A simple templates for MonteCarlo Simulation using MPI packages, using IsingModel as example.
# How to use it?
Inside the __main.c__,you will find  a function named: __monte_carlo_sampling()__ , which acts as a single-processor MonteCarlo program. When you change your model from 2D ising(in this example) to others(such as XY-model,Hubbard model), all you need to change is to copy your original code inside the sub function __monte_carlo_sampling()__ , then you can use MPI to run a parellel version of your code.
  
The final data will be shown in both console and the file __MCdata.dat__ , I've upload an example data file here for you to check the format. The different coloumns means datas comes from different processors, while the rows indicates different rows.
# Email me
If you found BUGs or you have some suggestion about the code, you are always welcome to email me at: __ienore@mail.ustc.edu.cn__ .
