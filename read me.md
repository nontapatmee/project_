2D polluted jet in larminar cross flow channal 
  In main function
   1.Adjust nx and ny for number of node cell in x-axis and y-axis conseqeunty.
   2.Re => Re-number , dt => duration per a time step 
   3.Multiply coefficient in dx term for increasing x-axis'grid size
   4.function visualize shows array of values
   5.Adjust n_max in "for(int n=1;n<=n_max;n++)" line to determine max time steps simulation.
   
  In simulation_p function
   1.it_max => max number of iteration
   2.eps => norm limit or absolute tolerance eps ,rit => discrete L^2-norm
  
  
3D 
In main function
   1.Adjust nx, ny and nz for number of node cell in x-axis, y-axis and z-axis conseqeunty.
   2.Re => Re-number , dt => duration per a time step 
   3.Multiply coefficient in dx term for increasing x-axis'grid size
   4.function visualize shows array of values
   5.Adjust n_max in "for(int n=1;n<=n_max;n++)" line to determine max time steps simulation.
   
  In simulation_p function
   1.it_max => max number of iteration
   2.eps => norm limit or absolute tolerance eps ,rit => discrete L^2-norm
