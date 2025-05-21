within n5GDHC.components;

model hn_tubes
  // Additional type definitions
  type HeatOfFusion = Real(final quantity = "Heat of fusion", final unit = "J/kg");
  type VolumeHeatCapacity = Real(final quantity = "Volume heat capacity", final unit = "J.m-3.K-1)");
  import Modelica.Units.Conversions.{from_deg,to_deg};
  // Discretization
  parameter Integer N_s = 2 "Number of energy balances for surrounding soil";
  // States and initials
  parameter Modelica.Units.NonSI.Temperature_degC T_f_init = 10.0 "Initial temperature of fluid";
  parameter Modelica.Units.NonSI.Temperature_degC T_t_init = T_f_init "Initial temperature of tube";
  parameter Modelica.Units.NonSI.Temperature_degC T_s_0 = 13.5 "Initial soil temperature at that time at nearest DWD station";
  parameter Modelica.Units.NonSI.Temperature_degC T_s_bd_min = 5.0 "Minimum temperature of model soil boundary throughout a year";
  parameter Modelica.Units.NonSI.Temperature_degC T_s_bd_max = 20.0 "Maximum temperature of model soil boundary throughout a year";
  parameter Modelica.Units.NonSI.Temperature_degC T_s_init[N_s] = linspace(T_s_0, T_s_bd_min + ((-0.5*cos(2*Modelica.Constants.pi*(t0/3600.0 - T_s_bd_shift)/8760.0)) + 0.5)*(T_s_bd_max - T_s_bd_min), N_s) "Initial temperature of soil";
  parameter Real T_s_bd_shift = 900.0 "Value to shift T_s_bd according to DWD measurements";
  Modelica.Units.NonSI.Temperature_degC T_f_sup(start = T_f_init, fixed = true) "Temperature of fluid supply line";
  Modelica.Units.NonSI.Temperature_degC T_f_ret(start = T_f_init, fixed = true) "Temperature of fluid return line";
  Modelica.Units.NonSI.Temperature_degC T_t_sup(start = T_t_init, fixed = true) "Temperature of tube supply line";
  Modelica.Units.NonSI.Temperature_degC T_t_ret(start = T_t_init, fixed = true) "Temperature of tube return line";
  Modelica.Units.NonSI.Temperature_degC T_s_sup_1[N_s](start = T_s_init, each fixed = true) "Temperatures of soil supply line; start = T_s_init";
  Modelica.Units.NonSI.Temperature_degC T_s_sup_2[N_s](start = T_s_init, each fixed = true) "Temperatures of soil supply line; start = T_s_init";
  Modelica.Units.NonSI.Temperature_degC T_s_ret_1[N_s](start = T_s_init, each fixed = true) "Temperatures of soil return line; start = T_s_init";
  Modelica.Units.NonSI.Temperature_degC T_s_ret_2[N_s](start = T_s_init, each fixed = true) "Temperatures of soil return line; start = T_s_init";
  parameter Modelica.Units.SI.Time t0 = 0.0 "Time of the year (in seconds) at which the simulation starts (required for soil temperature boundary condition)";
  // Material properties
  parameter Real ratio_s_w = 0.10 "Relative water content of soil";
  parameter HeatOfFusion dH_w = 333.55e3 "Heat of fusion water";
  parameter Modelica.Units.SI.Density rho_f = 1.0e3 "Density fluid";
  parameter Modelica.Units.SI.Density rho_t = 0.96e3 "Density tube material (polyethylen)";
  parameter Modelica.Units.SI.Density rho_s = 2.0e3 "Density soil";
  parameter Modelica.Units.SI.SpecificHeatCapacity c_f = 3.8e3 "Specific heat capactiy fluid";
  parameter Modelica.Units.SI.SpecificHeatCapacity c_t = 1.9e3 "Specific heat capactiy tube";
  parameter Modelica.Units.SI.SpecificHeatCapacity c_s_d = 0.8e3 "Specific heat capacity soil, dry";
  parameter Modelica.Units.SI.SpecificHeatCapacity c_w = 4.18e3 "Specific heat capactiy water";
  parameter Modelica.Units.SI.SpecificHeatCapacity c_ice = 2.06e3 "Specific heat capacity ice";
  parameter VolumeHeatCapacity c_s_f = rho_s*((1 - ratio_s_w)*c_s_d + c_w*ratio_s_w) "Specific heat capacity (volume) soil, fluid";
  parameter VolumeHeatCapacity c_s_m = c_s_f + rho_s*(ratio_s_w*dH_w) "Specific heat capactiy (volume) soil, water mixed";
  parameter VolumeHeatCapacity c_s_s = rho_s*((1 - ratio_s_w)*c_s_d + c_ice*ratio_s_w) "Specific heat capacity (volume) soil, solid (iced)";
  parameter Modelica.Units.NonSI.Temperature_degC T_w_fluid = 0.0 "Lower temperature bound water liquid";
  parameter Modelica.Units.NonSI.Temperature_degC T_w_solid = -1.0 "Upper temperature bound water solid";
  parameter Modelica.Units.SI.CoefficientOfHeatTransfer alpha_t = 2.5e3 "Heat transfer coeff. tube";
  parameter Modelica.Units.SI.ThermalConductivity lambda_t = 0.41 "Therm. conductivity tube";
  parameter Modelica.Units.SI.ThermalConductivity lambda_s = 1.0 "Therm. conductivity soil";
  parameter Modelica.Units.SI.ThermalConductivity lambda_s_i = 2.18 "Therm. conductivity ice";
  // Geometrical propteries
  parameter Modelica.Units.SI.Radius r_t = 0.07 "Inner radius of network tube";
  parameter Modelica.Units.SI.Thickness t_t = 0.005 "Thickness of the tube";
  parameter Modelica.Units.SI.Thickness t_s = 1.0 "Thickness of surrounding soil";
  parameter Modelica.Units.SI.Length l = 1.0 "Length of the tube";
  parameter Modelica.Units.SI.Length d_t = 0.5 "Distance between supply and return tube (Info Hr. Wickenheißer)";
  parameter Modelica.Units.SI.Radius r_b = d_t/2 "Treshold for soil layers";
  Modelica.Units.SI.Length d_sr_exch "Distance of considered return and supply exchange are";
  parameter Modelica.Units.SI.Length d_12_exch = 0.1 "Distance of exchange area between part 1 and 2 of supply and return side";
  // Calculated variables
  Modelica.Units.SI.Radius r_s[N_s] "Radius of the respective soil layers";
  Modelica.Units.SI.Length cl[N_s] "Circular chord length of soil circle segments";
  Modelica.Units.SI.Length h[N_s] "Sagitta of soil circle element";
  Modelica.Units.SI.Length s[N_s] "arc length";
  Modelica.Units.SI.Angle beta "Angle of outer circular segment";
  parameter Modelica.Units.SI.Angle beta_circle = from_deg(360);
  Real adj_1[N_s] "Adjustment factor for Qdot_s_s (depending on circumference of remaining soil layer)";
  Real adj_2[N_s] "Adjustment factor for Qdot_s_s (depending on circumference of remaining soil layer)";
  Modelica.Units.SI.Area A_f "Cross sectional area pipe";
  Modelica.Units.SI.Volume V_f "Fluid volume";
  Modelica.Units.SI.Mass m_f "Fluid mass";
  Modelica.Units.SI.Volume V_t "Tube volume";
  Modelica.Units.SI.Mass m_t "Tube mass";
  Modelica.Units.SI.Area A_t "Inner area tube";
  Modelica.Units.SI.Volume V_s_1[N_s] "Soil volumes part 1";
  Modelica.Units.SI.Volume V_s_2[N_s] "Soil volumes part 2";
  Modelica.Units.SI.Mass m_s_1[N_s] "Soil masses part 1";
  Modelica.Units.SI.Mass m_s_2[N_s] "Soil masses part 2";
  Modelica.Units.SI.Area A_s_1[N_s] "Soil area part 1";
  Modelica.Units.SI.Area A_s_2[N_s] "Soil area part 2";
  Modelica.Units.SI.Area A_s_hollow[N_s] "Soil area of hollow cylinders";
  Modelica.Units.SI.Area A_12_exch[N_s] "Soil area for exchange between part 1 and part 2";
  Modelica.Units.SI.Area A_sr_exch[N_s] "Soil area for exchange between supply and return side";
  Modelica.Units.SI.Length cir[N_s] "Circumference soil layers";
  VolumeHeatCapacity c_s_sup_1[N_s] "Specific heat capactiy soil supply line part 1";
  VolumeHeatCapacity c_s_ret_1[N_s] "Specific heat capactiy soil return line part 1";
  VolumeHeatCapacity c_s_sup_2[N_s] "Specific heat capactiy soil supply line part 2";
  VolumeHeatCapacity c_s_ret_2[N_s] "Specific heat capactiy soil return line part 2";
  Modelica.Units.NonSI.Temperature_degC T_s_bd "Temperature of soil boundary";
  Modelica.Units.SI.HeatFlowRate Qdot_t_f_sup "Heat exchange tube-fluid supply line";
  Modelica.Units.SI.HeatFlowRate Qdot_t_f_ret "Heat exchange tube-fluid return line";
  Modelica.Units.SI.HeatFlowRate Qdot_t_s_sup "Heat exchange tube-soil supply line";
  Modelica.Units.SI.HeatFlowRate Qdot_t_s_ret "Heat exchange tube-soil return line";
  Modelica.Units.SI.HeatFlowRate Qdot_s_s_sup_1[N_s] "Heat exchange soil-soil supply line part 1";
  Modelica.Units.SI.HeatFlowRate Qdot_s_s_ret_1[N_s] "Heat exchange soil-soil return line part 1";
  Modelica.Units.SI.HeatFlowRate Qdot_s_s_sup_2[N_s] "Heat exchange soil-soil supply line part 2";
  Modelica.Units.SI.HeatFlowRate Qdot_s_s_ret_2[N_s] "Heat exchange soil-soil return line part 2";
  Modelica.Units.SI.HeatFlowRate Qdot_sup_exch[N_s] "Heat exchange soil between part 1 and part 2 of supply side";
  Modelica.Units.SI.HeatFlowRate Qdot_ret_exch[N_s] "Heat exchange soil between part 1 and part 2 of return side";
  Modelica.Units.SI.HeatFlowRate Qdot_s_s_exch[N_s] "Heat exchange soil boundary between supply and return side";
  Modelica.Units.SI.Energy Q_t_s_sup(start = 0, fixed = true);
  Modelica.Units.SI.Energy Q_t_s_ret(start = 0, fixed = true);
  Modelica.Units.SI.Energy Q_t_f_sup(start = 0, fixed = true);
  Modelica.Units.SI.Energy Q_t_f_ret(start = 0, fixed = true);
  Modelica.Units.SI.Energy Q_s_s_sup_1[N_s](each start = 0, each fixed = true);
  Modelica.Units.SI.Energy Q_s_s_ret_1[N_s](each start = 0, each fixed = true);
  Modelica.Units.SI.Energy Q_s_s_sup_2[N_s](each start = 0, each fixed = true);
  Modelica.Units.SI.Energy Q_s_s_ret_2[N_s](each start = 0, each fixed = true);
  Modelica.Units.SI.Energy Q_sup_exch[N_s](each start = 0, each fixed = true);
  Modelica.Units.SI.Energy Q_ret_exch[N_s](each start = 0, each fixed = true);
  Modelica.Units.SI.Energy Q_s_s_exch[N_s](each start = 0, each fixed = true);
  Modelica.Units.SI.PressureDifference delta_p_sup "Pressure drop within discretized tube length";
  Modelica.Units.SI.PressureDifference delta_p_ret "Pressure drop within discretized tube length";
  Real f_d_sup "Darcy friction factor";
  Real f_d_ret "Darcy friction factor";
  Real Re_sup "Reynolds number";
  Real Re_ret "Reynolds number";
  Modelica.Units.SI.DynamicViscosity mu_sup "Dynamic viscosity of fluid inside tubes";
  Modelica.Units.SI.DynamicViscosity mu_ret "Dynamic viscosity of fluid inside tubes";
  parameter String path_mu_lookup_inp = Modelica.Utilities.Files.loadResource("modelica://n5GDHC/input/Lookup_table_mu.txt");
  connectors.fluid_p inlet_hn_supply annotation(
    Placement(visible = true, transformation(origin = {0, 0}, extent = {{58, -10}, {78, 10}}, rotation = 0), iconTransformation(origin = {-134, -78}, extent = {{64, 28}, {84, 48}}, rotation = 0)));
  connectors.fluid_p outlet_hn_supply annotation(
    Placement(visible = true, transformation(origin = {0, 0}, extent = {{-78, -10}, {-58, 10}}, rotation = 0), iconTransformation(origin = {120, -78}, extent = {{-70, 28}, {-50, 48}}, rotation = 0)));
  connectors.fluid_p outlet_hn_return annotation(
    Placement(visible = true, transformation(origin = {68, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  connectors.fluid_p inlet_hn_return annotation(
    Placement(visible = true, transformation(origin = {-68, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Tables.CombiTable1Ds mu_lookup_sup(tableOnFile = true, tableName = "mu", fileName = path_mu_lookup_inp, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, extrapolation = Modelica.Blocks.Types.Extrapolation.HoldLastPoint) annotation(
    Placement(transformation(origin = {-28, 44}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Tables.CombiTable1Ds mu_lookup_ret(extrapolation = Modelica.Blocks.Types.Extrapolation.HoldLastPoint, fileName = path_mu_lookup_inp, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, tableName = "mu", tableOnFile = true) annotation(
    Placement(transformation(origin = {10, 44}, extent = {{-10, -10}, {10, 10}})));
equation
// Geometry
  A_f = Modelica.Constants.pi*r_t^2;
  V_f = Modelica.Constants.pi*r_t^2*l;
  m_f = rho_f*V_f;
  A_t = 2*Modelica.Constants.pi*r_t*l;
  V_t = Modelica.Constants.pi*((r_t + t_t)^2 - r_t^2)*l;
  m_t = rho_t*V_t;
  for k in 1:N_s loop
    r_s[k] = r_t + t_t + k*(t_s/N_s);
    h[k] = max(r_s[k] - r_b, 0);
    cl[k] = 2*sqrt(2*r_s[k]*h[k] - h[k]^2);
    s[k] = 2*r_s[k]*asin(cl[k]/(2*r_s[k]));
    cir[k] = 2*Modelica.Constants.pi*r_s[k];
    if k == 1 then
      A_s_1[k] = A_s_hollow[k]*(beta_circle - beta)/beta_circle;
      A_s_2[k] = (Modelica.Constants.pi*(r_s[k]^2)*beta/beta_circle - (r_s[k]^2*asin(cl[k]/(2*r_s[k])) - (cl[k]*(r_s[k] - h[k])/2))) - Modelica.Constants.pi*((r_t + t_t)^2)*beta/beta_circle;
      A_s_hollow[k] = Modelica.Constants.pi*(r_s[k]^2) - Modelica.Constants.pi*((r_t + t_t)^2);
    else
      A_s_1[k] = A_s_hollow[k]*(beta_circle - beta)/beta_circle;
      A_s_2[k] = (Modelica.Constants.pi*(r_s[k]^2)*beta/beta_circle - (r_s[k]^2*asin(cl[k]/(2*r_s[k])) - (cl[k]*(r_s[k] - h[k])/2))) - (Modelica.Constants.pi*(r_s[k - 1]^2)*beta/beta_circle - (r_s[k - 1]^2*asin(cl[k - 1]/(2*r_s[k - 1])) - (cl[k - 1]*(r_s[k - 1] - h[k - 1])/2)));
      A_s_hollow[k] = Modelica.Constants.pi*(r_s[k]^2) - Modelica.Constants.pi*(r_s[k - 1]^2);
    end if;
    adj_1[k] = A_s_1[k]/A_s_hollow[k];
    adj_2[k] = A_s_2[k]/A_s_hollow[k];
    V_s_1[k] = A_s_1[k]*l;
    V_s_2[k] = A_s_2[k]*l;
    m_s_1[k] = (rho_s*V_s_1[k]);
    m_s_2[k] = (rho_s*V_s_2[k]);
  end for;
  beta = 2*asin(cl[N_s]/(2*r_s[N_s]));
// Heat flows
  Qdot_t_f_sup = 2*Modelica.Constants.pi*l / (1/(r_t*alpha_t) + 1/lambda_t * log((t_t/2+r_t)/r_t)) * (T_f_sup - T_t_sup);
  Qdot_t_f_ret = 2*Modelica.Constants.pi*l / (1/(r_t*alpha_t) + 1/lambda_t * log((t_t/2+r_t)/r_t)) * (T_f_ret - T_t_ret);
  Qdot_t_s_sup = 2*Modelica.Constants.pi*l / (1/lambda_t*log((r_t+t_t)/(r_t+t_t/2)) + 1/lambda_s*log((r_t+t_t+t_s/N_s/2)/(r_t+t_t))) * (T_t_sup - T_s_sup_1[1]*adj_1[1] - T_s_sup_2[1]*(1 - adj_1[1]));
  Qdot_t_s_ret = 2*Modelica.Constants.pi*l / (1/lambda_t*log((r_t+t_t)/(r_t+t_t/2)) + 1/lambda_s*log((r_t+t_t+t_s/N_s/2)/(r_t+t_t))) * (T_t_ret - T_s_ret_1[1]*adj_1[1] - T_s_ret_2[1]*(1 - adj_1[1]));
// Energy balances
  m_f*c_f*der(T_f_sup) = inlet_hn_supply.mdot*c_f*(inlet_hn_supply.T - T_f_sup) - Qdot_t_f_sup;
  m_f*c_f*der(T_f_ret) = inlet_hn_return.mdot*c_f*(inlet_hn_return.T - T_f_ret) - Qdot_t_f_ret;
  m_t*c_t*der(T_t_sup) = Qdot_t_f_sup - Qdot_t_s_sup;
  m_t*c_t*der(T_t_ret) = Qdot_t_f_ret - Qdot_t_s_ret;
// Section 1
  for k in 1:N_s - 1 loop
    if k <= 2 then
      Qdot_s_s_sup_1[k] = (2*lambda_s*Modelica.Constants.pi*l/log((r_t + t_t + (k+0.5)*(t_s/N_s))/(r_t + t_t + (k - 0.5)*(t_s/N_s)))*(T_s_sup_1[k] - T_s_sup_1[k + 1]))*adj_1[k];
      Qdot_s_s_ret_1[k] = (2*lambda_s*Modelica.Constants.pi*l/log((r_t + t_t + (k+0.5)*(t_s/N_s))/(r_t + t_t + (k - 0.5)*(t_s/N_s)))*(T_s_ret_1[k] - T_s_ret_1[k + 1]))*adj_1[k];
    else
      Qdot_s_s_sup_1[k] = (2*lambda_s*Modelica.Constants.pi*l/log((r_t + t_t + (k+0.5)*(t_s/N_s))/(r_t + t_t + (k - 0.5)*(t_s/N_s)))*(T_s_sup_1[k] - T_s_sup_1[k + 1]))*adj_1[k];
      Qdot_s_s_ret_1[k] = (2*lambda_s*Modelica.Constants.pi*l/log((r_t + t_t + (k+0.5)*(t_s/N_s))/(r_t + t_t + (k - 0.5)*(t_s/N_s)))*(T_s_ret_1[k] - T_s_ret_1[k + 1]))*adj_1[k];
    end if;
  end for;
  Qdot_s_s_sup_1[N_s] = (2*lambda_s*Modelica.Constants.pi*l/log((r_t + t_t + N_s*(t_s/N_s))/(r_t + t_t + (N_s - 0.5)*(t_s/N_s)))*(T_s_sup_1[N_s] - T_s_bd))*adj_1[N_s];
  Qdot_s_s_ret_1[N_s] = (2*lambda_s*Modelica.Constants.pi*l/log((r_t + t_t + N_s*(t_s/N_s))/(r_t + t_t + (N_s - 0.5)*(t_s/N_s)))*(T_s_ret_1[N_s] - T_s_bd))*adj_1[N_s];
  for k in 1:N_s loop
    if T_s_sup_1[k] >= T_w_fluid then
      c_s_sup_1[k] = c_s_f;
    elseif T_s_sup_1[k] <= T_w_solid then
      c_s_sup_1[k] = c_s_s;
    else
      c_s_sup_1[k] = c_s_m;
    end if;
  end for;
  for k in 1:N_s loop
    if T_s_ret_1[k] >= T_w_fluid then
      c_s_ret_1[k] = c_s_f;
    elseif T_s_ret_1[k] <= T_w_solid then
      c_s_ret_1[k] = c_s_s;
    else
      c_s_ret_1[k] = c_s_m;
    end if;
  end for;
  V_s_1[1]*c_s_sup_1[1]*der(T_s_sup_1[1]) = Qdot_t_s_sup*adj_1[1] - Qdot_s_s_sup_1[1] - Qdot_sup_exch[1];
  V_s_1[1]*c_s_ret_1[1]*der(T_s_ret_1[1]) = Qdot_t_s_ret*adj_1[1] - Qdot_s_s_ret_1[1] - Qdot_ret_exch[1];
  for k in 2:N_s - 1 loop
    V_s_1[k]*c_s_sup_1[k]*der(T_s_sup_1[k]) = Qdot_s_s_sup_1[k - 1] - Qdot_s_s_sup_1[k] - Qdot_sup_exch[k];
    V_s_1[k]*c_s_ret_1[k]*der(T_s_ret_1[k]) = Qdot_s_s_ret_1[k - 1] - Qdot_s_s_ret_1[k] - Qdot_ret_exch[k];
  end for;
  V_s_1[N_s]*c_s_sup_1[N_s]*der(T_s_sup_1[N_s]) = Qdot_s_s_sup_1[N_s - 1] - Qdot_s_s_sup_1[N_s] - Qdot_sup_exch[N_s];
  V_s_1[N_s]*c_s_ret_1[N_s]*der(T_s_ret_1[N_s]) = Qdot_s_s_ret_1[N_s - 1] - Qdot_s_s_ret_1[N_s] - Qdot_ret_exch[N_s];
// Section 2
  for k in 1:N_s - 1 loop
    if k <= 2 then
      Qdot_s_s_sup_2[k] = (2*lambda_s*Modelica.Constants.pi*l/log((r_t + t_t + (k+0.5)*(t_s/N_s))/(r_t + t_t + (k - 0.5)*(t_s/N_s)))*(T_s_sup_2[k] - T_s_sup_2[k + 1]))*adj_2[k];
      Qdot_s_s_ret_2[k] = (2*lambda_s*Modelica.Constants.pi*l/log((r_t + t_t + (k+0.5)*(t_s/N_s))/(r_t + t_t + (k - 0.5)*(t_s/N_s)))*(T_s_ret_2[k] - T_s_ret_2[k + 1]))*adj_2[k];
    else
      Qdot_s_s_sup_2[k] = (2*lambda_s*Modelica.Constants.pi*l/log((r_t + t_t + (k+0.5)*(t_s/N_s))/(r_t + t_t + (k - 0.5)*(t_s/N_s)))*(T_s_sup_2[k] - T_s_sup_2[k + 1]))*adj_2[k];
      Qdot_s_s_ret_2[k] = (2*lambda_s*Modelica.Constants.pi*l/log((r_t + t_t + (k+0.5)*(t_s/N_s))/(r_t + t_t + (k - 0.5)*(t_s/N_s)))*(T_s_ret_2[k] - T_s_ret_2[k + 1]))*adj_2[k];
    end if;
  end for;
  Qdot_s_s_sup_2[N_s] = 0.0;
  Qdot_s_s_ret_2[N_s] = 0.0;
  for k in 1:N_s loop
    if T_s_sup_2[k] >= T_w_fluid then
      c_s_sup_2[k] = c_s_f;
    elseif T_s_sup_2[k] <= T_w_solid then
      c_s_sup_2[k] = c_s_s;
    else
      c_s_sup_2[k] = c_s_m;
    end if;
  end for;
  for k in 1:N_s loop
    if T_s_ret_2[k] >= T_w_fluid then
      c_s_ret_2[k] = c_s_f;
    elseif T_s_ret_2[k] <= T_w_solid then
      c_s_ret_2[k] = c_s_s;
    else
      c_s_ret_2[k] = c_s_m;
    end if;
  end for;
  V_s_2[1]*c_s_sup_2[1]*der(T_s_sup_2[1]) = Qdot_t_s_sup*adj_2[1] - Qdot_s_s_sup_2[1] + Qdot_sup_exch[1] - Qdot_s_s_exch[1];
  V_s_2[1]*c_s_ret_2[1]*der(T_s_ret_2[1]) = Qdot_t_s_ret*adj_2[1] - Qdot_s_s_ret_2[1] + Qdot_ret_exch[1] + Qdot_s_s_exch[1];
  for k in 2:N_s - 1 loop
    V_s_2[k]*c_s_sup_2[k]*der(T_s_sup_2[k]) = Qdot_s_s_sup_2[k - 1] - Qdot_s_s_sup_2[k] + Qdot_sup_exch[k] - Qdot_s_s_exch[k];
    V_s_2[k]*c_s_ret_2[k]*der(T_s_ret_2[k]) = Qdot_s_s_ret_2[k - 1] - Qdot_s_s_ret_2[k] + Qdot_ret_exch[k] + Qdot_s_s_exch[k];
  end for;
  V_s_2[N_s]*c_s_sup_2[N_s]*der(T_s_sup_2[N_s]) = Qdot_s_s_sup_2[N_s - 1] + Qdot_sup_exch[N_s] - Qdot_s_s_exch[N_s];
  V_s_2[N_s]*c_s_ret_2[N_s]*der(T_s_ret_2[N_s]) = Qdot_s_s_ret_2[N_s - 1] + Qdot_ret_exch[N_s] + Qdot_s_s_exch[N_s];
//
  d_sr_exch = t_s/2*sin(beta/2);
  for k in 1:N_s loop
    if k == 1 then
      A_sr_exch[k] = cl[k]*l;
    else
      A_sr_exch[k] = (cl[k] - cl[k - 1])*l;
    end if;
    if h[k] == 0 then
      Qdot_s_s_exch[k] = 0.0;
    else
      Qdot_s_s_exch[k] = lambda_s*A_sr_exch[k]/d_sr_exch*(T_s_sup_2[k] - T_s_ret_2[k]);
    end if;
  end for;
//
  for k in 1:N_s loop
    A_12_exch[k] = t_s/N_s*l;
    Qdot_sup_exch[k] = lambda_s*A_12_exch[k]/d_12_exch*(T_s_sup_1[k] - T_s_sup_2[k]);
    Qdot_ret_exch[k] = lambda_s*A_12_exch[k]/d_12_exch*(T_s_ret_1[k] - T_s_ret_2[k]);
  end for;
  T_s_bd = T_s_bd_min + ((-0.5*cos(2*Modelica.Constants.pi*((t0 + time)/3600.0 - T_s_bd_shift)/8760.0)) + 0.5)*(T_s_bd_max - T_s_bd_min);
  outlet_hn_supply.T = T_f_sup;
  outlet_hn_return.T = T_f_ret;
  delta_p_sup = (8*l*inlet_hn_supply.mdot^2)/(rho_f*Modelica.Constants.pi^2*(r_t*2)^5)*f_d_sup;
  delta_p_ret = (8*l*inlet_hn_return.mdot^2)/(rho_f*Modelica.Constants.pi^2*(r_t*2)^5)*f_d_ret;
  f_d_sup = if noEvent(Re_sup > 0) then 0.3164/(Re_sup^0.25) else 0.0;
// event suppresion due to numerical issues
  f_d_ret = if noEvent(Re_ret > 0) then 0.3164/(Re_ret^0.25) else 0.0;
// event suppresion due to numerical issues
  Re_sup = (inlet_hn_supply.mdot*r_t*2)/(mu_sup*A_f);
  Re_ret = (inlet_hn_return.mdot*r_t*2)/(mu_ret*A_f);
  mu_lookup_sup.u = inlet_hn_supply.T;
  mu_sup = mu_lookup_sup.y[1];
  mu_lookup_ret.u = inlet_hn_return.T;
  mu_ret = mu_lookup_ret.y[1];
  outlet_hn_supply.p = inlet_hn_supply.p - delta_p_sup;
  outlet_hn_return.p = inlet_hn_return.p - delta_p_ret;
// Mass balance
  outlet_hn_supply.mdot = inlet_hn_supply.mdot;
  outlet_hn_return.mdot = inlet_hn_return.mdot;
// Monitoring
  der(Q_s_s_sup_1) = Qdot_s_s_sup_1;
  der(Q_s_s_ret_1) = Qdot_s_s_ret_1;
  der(Q_s_s_sup_2) = Qdot_s_s_sup_2;
  der(Q_s_s_ret_2) = Qdot_s_s_ret_2;
  der(Q_t_f_sup) = Qdot_t_f_sup;
  der(Q_t_f_ret) = Qdot_t_f_ret;
  der(Q_t_s_sup) = Qdot_t_s_sup;
  der(Q_t_s_ret) = Qdot_t_s_ret;
  der(Q_s_s_exch) = Qdot_s_s_exch;
  der(Q_sup_exch) = Qdot_sup_exch;
  der(Q_ret_exch) = Qdot_ret_exch;
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-90, -80}, {90, 80}}), graphics = {Rectangle(origin = {8, 40}, lineThickness = 0.5, extent = {{-68, 30}, {52, -30}}), Ellipse(origin = {8, 40}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{34, 30}, {68, -30}}), Line(origin = {6.78479, 39.8101}, points = {{-36, 0}, {26, 0}}, thickness = 0.5), Line(origin = {6.78479, 39.8101}, points = {{-42, 0}, {-36, 6}, {-36, -6}, {-42, 0}}, thickness = 0.5), Rectangle(origin = {8, -40}, lineThickness = 0.5, extent = {{-68, 30}, {52, -30}}), Ellipse(origin = {-110, -40}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{34, 30}, {68, -30}}), Line(origin = {2.71745, -40.0193}, points = {{-36, 0}, {26, 0}}, thickness = 0.5), Line(origin = {65.2707, -40.0193}, points = {{-30, 0}, {-36, 6}, {-36, -6}, {-30, 0}}, thickness = 0.5), Ellipse(origin = {-111, 40}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{34, 30}, {68, -30}}, startAngle = 90, endAngle = 270, closure = EllipseClosure.Chord), Ellipse(origin = {-9, 40}, rotation = 180, fillColor = {255, 255, 255}, pattern = LinePattern.Dash, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{34, 30}, {68, -30}}, startAngle = 90, endAngle = 270, closure = EllipseClosure.Chord), Ellipse(origin = {-22, 40}, rotation = 180, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{25, 29.75}, {49, -29.75}}), Ellipse(origin = {9, -40}, fillColor = {255, 255, 255}, pattern = LinePattern.Dash, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{34, 30}, {68, -30}}, startAngle = 90, endAngle = 270, closure = EllipseClosure.Chord), Ellipse(origin = {111, -40}, rotation = 180, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{34, 30}, {68, -30}}, startAngle = 90, endAngle = 270, closure = EllipseClosure.Chord), Ellipse(origin = {96, -40}, rotation = 180, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{25, 29.75}, {49, -29.75}})}),
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-90, -80}, {90, 80}})));
end hn_tubes;
