function MORFE_integrate_rdyn_frc(analysis_name,
  zero_amplitude,harmonics_init,param_init,cont_param,time_integration_length,
  forward=true,MaxNumPoints=100,minstep=1e-8,maxstep=20.0,initstep=0.1,ncol=4.0,ntst=40.0,analysis_number=1)
    # 
    # check for continuation direction
    if (forward)
      var = 0
    else
      var = 1
    end
    # 
    # output folder
    analysis_folder = "./"*analysis_name*"/matcont_automatic"
    analysis_output = "./"*analysis_name*"/frc_"*string(analysis_number)*".mat"
    # 
    nΩ = Int(size(harmonics_init)[1]/2)
    control_parameters = "mu"
    for i = 1:nΩ
      control_parameters *= ",beta"*string(i)
    end
    for i = 1:nΩ
      control_parameters *= ",w"*string(i)
    end
    # 
    # initialise normal coordinates and harmonics amplitudes
    X0 = zeros(Float64,size(zero_amplitude)[1]+size(harmonics_init)[1])
    X0[1:size(zero_amplitude)[1]] = zero_amplitude
    X0[size(zero_amplitude)[1]+1:end] = harmonics_init
    # 
    # inizialise control parameters
    param = ""
    param *= "mu = "*string(param_init[1])
    for i = 1:nΩ
      param *= ";beta"*string(i) * " = " * string(param_init[i+1])
    end
    for i = 1:nΩ
      param *= ";w"*string(i) * " = " * string(param_init[i+1+nΩ])
    end
    param *= ";"
    #
    # Matlab session
    ms = MSession();
    @mput analysis_folder
    @mput analysis_output
    @mput X0
    @mput time_integration_length
    @mput cont_param
    @mput ncol
    @mput ntst
    @mput MaxNumPoints
    @mput maxstep
    @mput minstep
    @mput initstep
    @mput var
    #
    eval_string(ms,"
    warning off;
    addpath(genpath('./MatCont7p3')) 
    addpath(analysis_folder);
    ndofs = size(X0);
    ndofs = ndofs(1);  
    init 
    tfin = time_integration_length;
    ")
    #
    eval_string(ms,param)
    #
    eval_string(ms,"
    %%           Find periodic limit cycle for first continuation point
    ndofs = size(X0);
    ndofs = ndofs(1);    
    hls = feval(@MORFEsystem);
    options=odeset('RelTol',1e-8);
    [t,y] = ode45(hls{2},[0,tfin],X0,options,"*control_parameters*");
    x1 = y(end,:);
    [vl,pk] = findpeaks(y(:,1));
    rel_err = (vl(end)-vl(end-1))/vl(end);
    period = t(pk(end))-t(pk(end-1));
    %
    %
    disp(\"-----------------------------------\")
    disp(\"| Time integration to find periodic limit cycle at first continuation point\")
    disp(\"| Relative error between last two consecutive peaks:\")    
    disp(strcat(\"| \",num2str(rel_err)))
    disp(\"| Frequency (rad,s) found at supposed steady state:\")    
    disp(strcat(\"| \",num2str(2*pi/period)))
    disp(\"| Consider enlarging the time_integration_length if the transient doesn't appear decayed\")
    disp(\"| Figure 100 shows the behaviour of the first normal coordinates in time\")
    figure(100);plot(t,y(:,1));
    xlabel('\$t\$','Interpreter','latex');
    ylabel('\$a_1 \$','Interpreter','latex');
    disp(\"-----------------------------------\")
    %
    %              integrate for a single period
    [t,y] = ode45(hls{2},[0 period],x1,options,"*control_parameters*");
    tolerance=1e-4;
    %
    %              use limit cycle as initial orbit for continuation
    [x0,v0]=initOrbLC(@MORFEsystem,t,y,["*control_parameters*"],[cont_param],ntst,ncol,tolerance);
    %
    %%           Continuation of Periodic orbits
    opt=contset; 
    opt=contset(opt,'MaxNumPoints',MaxNumPoints); 
    opt=contset(opt, 'InitStepsize' , initstep); 
    opt=contset(opt,'MaxStepsize'  , maxstep); 
    opt=contset(opt,'MinStepsize'  , minstep); 
    opt=contset(opt,'Backward',var); 
    opt=contset(opt,'FunTolerance', 1e-6);
    opt=contset(opt,'VarTolerance', 1e-6);
    opt=contset(opt,'ActiveParams', cont_param);
    %
    %              launch continuation
    [xlcc,vlcc,slcc,hlcc,flcc]=cont(@limitcycle,x0,v0,opt); 
    %
    %              save output
    save(analysis_output,'xlcc','ncol','ntst')
    ")
    #
    return ms
    #
  end


function MORFE_integrate_rdyn_backbone(analysis_name,
  X0,time_integration_length,forward=true,MaxNumPoints=100,
  minstep=1e-8,maxstep=20.0,initstep=0.1,ncol=4.0,ntst=40.0,analysis_number=1)

  if (forward)
    var = 0
  else
    var = 1
  end

  analysis_folder = "./"*analysis_name*"/matcont_automatic"
  analysis_output = "./"*analysis_name*"/backbone_"*string(analysis_number)*".mat"

  ms = MSession();
  @mput analysis_folder
  @mput analysis_output
  @mput X0
  @mput time_integration_length
  @mput ncol
  @mput ntst
  @mput MaxNumPoints
  @mput maxstep
  @mput minstep
  @mput initstep
  @mput var
  eval_string(ms,"
  warning off
  addpath(genpath('./MatCont7p3')) 
  addpath(analysis_folder)
  init 
  mu = 0.0;
  ndofs = size(X0);
  ndofs = ndofs(1);  
  tfin=time_integration_length; 
  %
  hls=feval(@MORFEsystem);
  options=odeset('RelTol',1e-8);
  %
  [t,y]=ode45(hls{2},[0,tfin],X0,options,mu);
  %
  x1 = y(end,:);
  %
  [vl,pk]=findpeaks(y(:,1));
  rel_err = (vl(end)-vl(end-1))/vl(end);
  period = t(pk(end))-t(pk(end-1));
  %
  %
  disp(\"-----------------------------------\")
  disp(\"| Time integration to find periodic limit cycle at first continuation point\")
  disp(\"| Relative error between last two consecutive peaks:\")    
  disp(strcat(\"| \",num2str(rel_err)))
  disp(\"| Frequency (rad,s) found at supposed steady state:\")    
  disp(strcat(\"| \",num2str(2*pi/period)))
  disp(\"| Consider enlarging the time_integration_length if the transient doesn't appear decayed\")
  disp(\"| Figure 100 shows the behaviour of the first normal coordinates in time\")
  figure(100);plot(t,y(:,1));
  xlabel('\$t\$','Interpreter','latex');
  ylabel('\$a_1 \$','Interpreter','latex');
  disp(\"-----------------------------------\")
  %
  %
  [t,y] = ode45(hls{2},[0 period],x1,options,mu);
  %
  %
  active_pars=[1]; 
  tolerance=1e-4;
  [x0,v0]=initOrbLC(@MORFEsystem,...
  t,y,...
  [mu],active_pars,ntst,ncol,...
  tolerance);
  %
  opt=contset; 
  opt=contset(opt,'MaxNumPoints',MaxNumPoints); 
  opt=contset(opt, 'InitStepsize' , initstep); 
  opt=contset(opt,'MaxStepsize'  , maxstep); 
  opt=contset(opt,'MinStepsize'  , minstep); 
  opt=contset(opt,'Backward',var); 
  opt=contset(opt,'FunTolerance', 1e-6);
  opt=contset(opt,'VarTolerance', 1e-6);
  [xlcc,vlcc,slcc,hlcc,flcc]=cont(@limitcycle,x0,v0,opt); 
  %
  save(analysis_output,'xlcc','ncol','ntst')
  ")
  #
  return ms
  #

end



function MORFE_compute_backbone_modal(analysis_name,analysis_number=1)
  #
  comb = readdlm(analysis_name*"/comb_a_0/monomials.txt")
  maps = readdlm(analysis_name*"/manifold_a_0/psi.txt")
  #
  po_dir = "./"*analysis_name*"/backbone_"*string(analysis_number)*".mat"
  #
  nc = size(comb)[1]
  ndofs = size(comb)[2]
  neig = size(maps)[2]
  #
  data = matread(po_dir)
  #
  xlcc = data["xlcc"]
  ncol = data["ncol"]
  ntst = data["ntst"]
  #
  npoints = Int(ncol*ntst)
  npo = size(xlcc)[2]
  vars = zeros(Float64,ndofs)
  field = zeros(Float64,neig)
  frf = zeros(npo,neig+1)
  #
  for po = 1:npo
  #
    ω = 2*π/xlcc[end-1,po]
    frf[po,1] = ω
    for p = 1:npoints
      # extract variables
      fill!(field,0.0)
      #
      for d = 1:ndofs
        vars[d] = xlcc[d+(p-1)*ndofs,po]
      end
      #
      for c = 1:nc
        #
        Πi = 1.0
        for d = 1:ndofs 
          Πi *= vars[d]^comb[c,d]
        end
        #
        for d = 1:neig
          field[d] += maps[c,d]*Πi
        end
        #
      end
      #
      for d = 1:neig
        if (abs(field[d])>frf[po,1+d])
          frf[po,1+d] = abs(field[d])
        end
      end
      #
    end
    #
  end
  #
  return frf
  #
end

# reconstruction with only autonomous vector till the paper on the forcing is
# becomes public
function MORFE_compute_frc_modal(analysis_name,Ω_list,analysis_number=1)
  #
  comb = readdlm(analysis_name*"/comb_a_0/monomials.txt")
  maps = readdlm(analysis_name*"/manifold_a_0/psi.txt")
  #
  po_dir = "./"*analysis_name*"/frc_"*string(analysis_number)*".mat"
  #
  nc = size(comb)[1]
  ndofs = size(comb)[2]
  nΩ = 2*size(Ω_list)[1]
  neig = size(maps)[2]
  #
  data = matread(po_dir)
  #
  xlcc = data["xlcc"]
  ncol = data["ncol"]
  ntst = data["ntst"]
  #
  npoints = Int(ncol*ntst)
  npo = size(xlcc)[2]
  vars = zeros(Float64,ndofs)
  field = zeros(Float64,neig)
  frf = zeros(npo,neig+1)
  #
  for po = 1:npo
  #
    ω = xlcc[end,po]
    frf[po,1] = ω
    for p = 1:npoints
      # extract variables
      fill!(field,0.0)
      #
      for d = 1:ndofs
        vars[d] = xlcc[d+(p-1)*(ndofs+nΩ),po]
      end
      #
      for c = 1:nc
        #
        Πi = 1.0
        for d = 1:ndofs 
          Πi *= vars[d]^comb[c,d]
        end
        #
        for d = 1:neig
          field[d] += maps[c,d]*Πi
        end
        #
      end
      #
      for d = 1:neig
        if (abs(field[d])>frf[po,1+d])
          frf[po,1+d] = abs(field[d])
        end
      end
      #
    end
    #
  end
  #
  return frf
  #
end