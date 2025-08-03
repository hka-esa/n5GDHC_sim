/*

This file is part of n5GDHC_sim.

Copyright (c) 2025 Manuel Kollmar, Adrian Bürger, Markus Bohlayer, Marco Braun, Moritz Diehl.
Developed at HS Karlsruhe and IMTEK, University of Freiburg.
All rights reserved.

The BSD 3-Clause License

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

within n5GDHC.components;

/* 
The parameter data is taken from literature, public sources or data sheets in order to reflect the specific conditions at the site.
Sources:
[1] Kuchling, H.: Taschenbuch der Physik, Hanser München, 22, 2022.
[2] Dehner, U.: Boden und Energiewende, Springer Fachmedien Wiesbaden, 2015.
[3] Herr, H.: Technische Physik: Wärmelehre, Europa-Lehrmittel, 2006.
[4] Deutscher Wetterdienst, Open Data, https://opendata.dwd.de/, 2024.
[5] Schweizer A.: Stoffdaten von Wasser in Abhängigkeit der Temperatur, https://schweizer-fn.de/stoff/wasser/wasser_stoff.php, 2025.
*/

model is
  // Functions and additional type definitions
  import Modelica.Units.Conversions.from_kWh;
  import Modelica.Units.Conversions.to_kWh;
  type HeatOfFusion = Real(final quantity = "Heat of fusion", final unit = "J/kg");
  type VolumeHeatCapacity = Real(final quantity = "Volume heat capacity", final unit = "J.m-3.K-1)");
  // Discretization
  parameter Integer N_hx = 30 "Number of heat exchanger strands";
  parameter Integer N_s = 10 "Number of energy balances for surrounding soil";
  parameter Integer N_w = 10 "Number of energy balances for water layers";
  // States and initials
  parameter Modelica.Units.NonSI.Temperature_degC T_w_1_init = 10.0 "Initial temperature water bottom layer";
  parameter Modelica.Units.NonSI.Temperature_degC T_w_Nw_init = 10.0 "Initial temperature water top layer";
  parameter Modelica.Units.NonSI.Temperature_degC T_w_init[N_w] = linspace(T_w_1_init, T_w_Nw_init, N_w);
  parameter Modelica.Units.NonSI.Temperature_degC T_hx_init = 10.0 "Initial heat exchanger temperature";
  parameter Modelica.Units.NonSI.Temperature_degC T_c_init[N_w] = T_w_init "Initial temperature concrete";
  parameter Modelica.Units.NonSI.Temperature_degC T_s_0 = 13.5 "Initial soil temperature at that time at nearest DWD station";
  parameter Modelica.Units.NonSI.Temperature_degC T_s_init[N_s] = linspace(T_s_0, T_s_bd_min + ((-0.5*cos(2*Modelica.Constants.pi*(t0/3600.0 - T_s_bd_shift)/8760.0)) + 0.5)*(T_s_bd_max - T_s_bd_min), N_s);
  parameter Real T_s_bd_shift = 900.0 "Value to shift T_s_bd according to DWD measurements";
  Modelica.Units.NonSI.Temperature_degC T_w[N_w](start = T_w_init, each fixed = true) "Temperature water";
  Modelica.Units.NonSI.Temperature_degC T_c[N_w](start = T_c_init, each fixed = true) "Temperature concrete";
  Modelica.Units.NonSI.Temperature_degC[N_w, N_s] T_s(each start = T_s_init[1], each fixed = false) "Temperature soil";
  Modelica.Units.NonSI.Temperature_degC T_hx[N_w](each start = T_hx_init, each fixed = true) "Temperature heat exchanger";
  // Boundary conditions
  parameter Modelica.Units.SI.Time t0 = 0.0 "Initial time of the year (in seconds) at which the simulation starts";
  parameter Modelica.Units.SI.Time t0_h = 28512000 "Start of heating period in seconds since simulation start";
  parameter Modelica.Units.SI.Time tN_h = 15292800 "End of heating period in seconds since simulation start";
  parameter Modelica.Units.NonSI.Temperature_degC T_s_bd_min = 5.0 "Minimum temperature of model soil boundary throughout a year [4]";
  parameter Modelica.Units.NonSI.Temperature_degC T_s_bd_max = 20.0 "Maximum temperature of model soil boundary throughout a year [4]";
  // Material properties
  parameter HeatOfFusion dH_w = 333.55e3 "Heat of fusion water [1]";
  constant Modelica.Units.NonSI.Temperature_degC T_w_l = 0.0 "Temperature above which all water is liquid";
  constant Modelica.Units.NonSI.Temperature_degC T_w_f = -1.0 "Temperature below which all water is frozen";
  parameter Modelica.Units.SI.Density rho_w = 1.0e3 "Density water [1]";
  parameter Modelica.Units.SI.Density rho_hx = 1.0e3 "Density brine [1]";
  parameter Modelica.Units.SI.Density rho_c = 2.0e3 "Density concrete [1]";
  parameter Modelica.Units.SI.Density rho_s = 2.0e3 "Density soil [1]";
  parameter Modelica.Units.SI.SpecificHeatCapacity c_w = 4.18e3 "Specific heat capactiy water [1]";
  parameter Modelica.Units.SI.SpecificHeatCapacity c_hx = 3.8e3 "Specific heat capactiy heat exchanger medium";
  parameter Modelica.Units.SI.SpecificHeatCapacity c_ice = 2.06e3 "Specific heat capacity ice [1]";
  parameter Modelica.Units.SI.SpecificHeatCapacity c_c = 0.84e3 "Specific heat capactiy concrete [1]";
  parameter Real ratio_s_w = 0.1 "Share of volumetric water in soil";
  parameter Modelica.Units.SI.SpecificHeatCapacity c_s_d = 0.8e3 "Specific heat capactiy soil dry [1,2]";
  parameter VolumeHeatCapacity c_s_f = ((1-ratio_s_w) * c_s_d + c_w * ratio_s_w) "Volumetric heat capacity soil with fluid water";
  parameter VolumeHeatCapacity c_s_m = c_s_f + rho_s * (ratio_s_w * dH_w) "Volumetric heat capacity soil while freezing";
  parameter VolumeHeatCapacity c_s_s = rho_s * ((1-ratio_s_w) * c_s_d + c_ice * ratio_s_w) "Volumetric heat capacity soil with frozen water";
  parameter Modelica.Units.SI.CoefficientOfHeatTransfer alpha_w_c = 400.0 "Heat transfer coefficient water - concrete [3]";
  parameter Modelica.Units.SI.CoefficientOfHeatTransfer alpha_b_hx = 4000.0 "Heat transfer coefficient brine - hx [1]";
  parameter Modelica.Units.SI.CoefficientOfHeatTransfer alpha_hx_w = 250.0 "Heat transfer coefficient hx - water [1]";
  parameter Modelica.Units.SI.ThermalConductivity lambda_c = 2.1 "Therm. conductivity concrete [1]";
  parameter Modelica.Units.SI.ThermalConductivity lambda_ice = 2.25 "Therm. conductivity ice [1]";
  parameter Modelica.Units.SI.ThermalConductivity lambda_t = 0.41 "Therm. conductivity tube [1]";
  parameter Modelica.Units.SI.ThermalConductivity lambda_s = 1.6 "Therm. conductivity soil [1]";
  parameter Modelica.Units.SI.ThermalConductivity lambda_w = 0.582 "Thermal Conductivity of water at around 10°C [5]";
  VolumeHeatCapacity c_s[N_w, N_s] "Volumetric heat capacity soil layers";
  // Geometrical properties
  parameter Modelica.Units.SI.Volume V_w = 543.0 "Water volume of the ice storage";
  parameter Modelica.Units.SI.Thickness t_c = 0.3 "Thickness of the concrete wall";
  parameter Modelica.Units.SI.Thickness t_s = 1.0 "Thickness of considered surrounding soil";
  parameter Modelica.Units.SI.Diameter d_t = 0.04 "Inner diameter of the tube";
  parameter Modelica.Units.SI.Thickness t_t = 0.004 "Thickness of the tube";
  parameter Modelica.Units.SI.Length l_t = 1e4 "HX length";
  parameter Modelica.Units.SI.Diameter d_w = 10.0 "Diameter of the water cylinder";
  Modelica.Units.SI.Volume V "Total volume of storage";
  Modelica.Units.SI.Mass m_w "Mass of a water volume";
  Modelica.Units.SI.Radius r_w "Radius of the cylinder holding the water volume";
  Modelica.Units.SI.Height h_w "Height of the cylinder holding the water volume";
  Modelica.Units.SI.Area A_w "Surface area of a water cylinder layer";
  Modelica.Units.SI.Volume V_hx "Heat exchanger volume";
  Modelica.Units.SI.Mass m_hx "Mass of a heat exchanger volume";
  Modelica.Units.SI.Area A_t "Cross-sectional area of heat exchanger tubes";
  Modelica.Units.SI.Area A_hx "Surface area of a heat exchanger volume";
  Modelica.Units.SI.MassFlowRate mdot_hx "Mass flow per heat exchanger volume";
  Modelica.Units.SI.Radius r_ice[N_w] "Radius of ice around heat exchanger";
  Modelica.Units.SI.Volume V_c[N_w] "Concrete volume";
  Modelica.Units.SI.Mass m_c[N_w] "Concrete mass";
  Modelica.Units.SI.Volume V_s[N_w, N_s] "Soil volume per discrete volume";
  Modelica.Units.SI.HeatFlowRate Qdot_hx[N_w] "Heat exchanger heat flow";
  Modelica.Units.SI.HeatFlowRate Qdot_m[N_w-1] "Heat exchange between storage volumes";
  Modelica.Units.SI.HeatFlowRate Qdot_w_c[N_w] "Heat flow concrete - water";
  Modelica.Units.SI.HeatFlowRate Qdot_c_s[N_w] "Heat flow concrete - soil";
  Modelica.Units.SI.HeatFlowRate Qdot_s_s[N_w, N_s - 1] "Heat flow soil - soil";
  Modelica.Units.SI.HeatFlowRate Qdot_s_bd[N_w] "Heat flow soil - soil boundary";
  Modelica.Units.SI.HeatCapacity C_w[N_w] "Heat capacity of water mass, depending on temperature";
  Modelica.Units.NonSI.Temperature_degC T_s_bd "Temperature of soil boundary";
  Real p_ice[N_w] "Degree of icing";
  Real x_ice[N_w](each start=0, each fixed=true) "Share of ice as state";
  Modelica.Units.SI.PressureDifference delta_p "Pressure drop within heat exchanger tubes";
  Real f_d "Darcy friction factor";
  Real Re "Reynolds number";
  Modelica.Units.SI.DynamicViscosity mu "Dynamic viscosity of fluid inside tubes";
  parameter String path_mu_lookup_inp = Modelica.Utilities.Files.loadResource("modelica://n5GDHC/input/Lookup_table_mu.txt");
  // Monitoring
  Modelica.Units.SI.Heat Q_hx(start = 0.0, fixed = true);
  Modelica.Units.SI.Heat Q_s(start = 0.0, fixed = true);
  Modelica.Units.SI.Heat Q_s_bd(start = 0.0, fixed = true);
  Modelica.Units.SI.HeatFlowRate Qdot_hx_sum;
  Modelica.Units.SI.HeatFlowRate Qdot_s_s_sum;
  Modelica.Units.SI.HeatFlowRate Qdot_c_s_sum;
  Modelica.Units.NonSI.Energy_kWh Q_hx_kWh;
  Modelica.Units.NonSI.Energy_kWh Q_s_kWh;
  connectors.fluid_p outlet_hx annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}}), iconTransformation(extent = {{50, -72}, {70, -52}})));
  connectors.fluid_p inlet_hx annotation(
    Placement(visible = true, iconTransformation(origin = {0, 0}, extent = {{50, 12}, {70, 32}}, rotation = 0)));
  Modelica.Blocks.Tables.CombiTable1Ds mu_lookup(tableOnFile = true, tableName = "mu", fileName = path_mu_lookup_inp, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, extrapolation = Modelica.Blocks.Types.Extrapolation.LastTwoPoints)  annotation(
    Placement(transformation(origin = {-32, 50}, extent = {{-10, -10}, {10, 10}})));
initial equation
  for k in 1:N_w loop
    T_s[k, :] = T_s_init;
  end for;
equation
// Geometry
  A_t = Modelica.Constants.pi*(d_t/2)^2;
  A_w = Modelica.Constants.pi*r_w^2;
  V_hx = Modelica.Constants.pi*(d_t/2)^2*l_t;
  V = V_w + V_hx;
  r_w = d_w/2;
  h_w = V/(Modelica.Constants.pi*r_w^2);
  m_w = rho_w*V_w/N_w;
  m_hx = V_hx*rho_hx/N_hx;
  A_hx = Modelica.Constants.pi*d_t*l_t/N_hx;
  mdot_hx = inlet_hx.mdot/N_hx;
  for k in 1:N_w loop
    if k == 1 or k == N_w then
      V_c[k] = Modelica.Constants.pi*((r_w + t_c)^2 - r_w^2)*(h_w/N_w) + Modelica.Constants.pi*r_w^2*t_c;
    else
      V_c[k] = Modelica.Constants.pi*((r_w + t_c)^2 - r_w^2)*(h_w/N_w);
    end if;
    m_c[k] = V_c[k]*rho_c;
    for l in 1:N_s loop
      if k == 1 or k == N_w then
        V_s[k, l] = Modelica.Constants.pi*((r_w + t_c + l*(t_s/N_s))^2 - (r_w + t_c + (l - 1)*(t_s/N_s))^2)*(h_w/N_w) + Modelica.Constants.pi*r_w^2*(t_s/N_s);
      else
        V_s[k, l] = Modelica.Constants.pi*((r_w + t_c + l*(t_s/N_s))^2 - (r_w + t_c + (l - 1)*(t_s/N_s))^2)*(h_w/N_w);
      end if;
    end for;
  end for;
  for k in 1:N_w loop
    for l in 1:N_s loop
      if T_s[k, l] >= T_w_l then
        c_s[k, l] = c_s_f;
      elseif T_s[k, l] <= T_w_f then
        c_s[k, l] = c_s_s;
      else
        c_s[k, l] = c_s_m;
      end if;
    end for;
  end for;
// Heat flows
  for k in 1:N_w-1 loop
    Qdot_m[k] = lambda_w * A_w / h_w * (T_w[k+1] - T_w[k]);
  end for;
  for k in 1:N_w loop
    if k == 1 or k == N_w then
      Qdot_w_c[k] = Modelica.Constants.pi / (1/((2*r_w*h_w/N_w + r_w^2)*alpha_w_c) + 1/(2*lambda_c*h_w/N_w)*log((r_w+t_c/2)/r_w) + t_c/2/(lambda_c*r_w^2)) * (T_w[k] - T_c[k]);
      Qdot_c_s[k] = Modelica.Constants.pi / (1/(2*lambda_c*h_w/N_w)*log((r_w+t_c)/(r_w+t_c/2)) + t_c/2/(lambda_c*r_w^2) + 1/(2*lambda_s*h_w/N_w)*log((r_w+t_c+t_s/N_s/2)/(r_w+t_c)) + (t_c+t_s/N_s/2)/(lambda_s*r_w^2)) * (T_c[k] - T_s[k,1]);
    else
      Qdot_w_c[k] = 2*Modelica.Constants.pi*(h_w/N_w) / ( 1/(r_w*alpha_w_c) + log((r_w + t_c/2)/(r_w))*1/lambda_c) * (T_w[k] - T_c[k]);
      Qdot_c_s[k] = 2*Modelica.Constants.pi*(h_w/N_w) / (1/lambda_c*log((r_w+t_c)/(r_w+t_c/2)) + 1/lambda_s*log((r_w+t_c+t_s/N_s/2)/(r_w+t_c))) * (T_c[k] - T_s[k,1]);
    end if;
    for l in 1:N_s - 1 loop
      if k == 1 or k == N_w then
        Qdot_s_s[k, l] = (2*lambda_s*Modelica.Constants.pi*(h_w/N_w)/log((r_w + t_c + (l+0.5)*(t_s/N_s))/(r_w + t_c + (l - 0.5)*(t_s/N_s))) + (lambda_s/(t_s/N_s))*(Modelica.Constants.pi*r_w^2))*(T_s[k, l] - T_s[k, l + 1]);
      else
        Qdot_s_s[k, l] = 2*lambda_s*Modelica.Constants.pi*(h_w/N_w)/log((r_w + t_c + (l+0.5)*(t_s/N_s))/(r_w + t_c + (l - 0.5)*(t_s/N_s)))*(T_s[k, l] - T_s[k, l + 1]);
      end if;
    end for;
    Qdot_s_bd[k] = 2*lambda_s*Modelica.Constants.pi*(h_w/N_w)/log((r_w + t_c + N_s*(t_s/N_s))/(r_w + t_c + (N_s - 0.5)*(t_s/N_s)))*(T_s[k, N_s] - T_s_bd);
// Energy balances
    assert(T_w[k] >= T_w_f, "Storage layer completely frozen", level = AssertionLevel.warning);
    if T_w[k] >= T_w_l then
      C_w[k] = m_w*c_w;
    elseif T_w[k] < T_w_l and T_w[k] >= T_w_f then
      C_w[k] = m_w*(dH_w*1.0 + c_w);
    else
      C_w[k] = m_w*c_ice;
    end if;
    r_ice[k] = sqrt((V_w/N_w * p_ice[k] + Modelica.Constants.pi * ((d_t/2)+t_t)^2 * l_t/N_hx) / (l_t / N_hx * Modelica.Constants.pi));
    if time < tN_h or time > t0_h then
      Qdot_hx[k] = 2 * Modelica.Constants.pi * (l_t/N_hx/N_w) / ( (1/(alpha_hx_w * (d_t/2))) + (log(r_ice[k] / ((d_t/2)+t_t)) / lambda_ice) + (log(((d_t/2)+t_t)/(d_t/2)) / lambda_t) + 1/(alpha_b_hx * (d_t/2)) ) * (T_hx[k] - T_w[k]);
    else
      Qdot_hx[k] = 2 * Modelica.Constants.pi * (l_t/N_hx/N_w) / ( (1/(alpha_hx_w * (d_t)/2)) + (log(((d_t/2)+t_t)/(d_t/2)) / lambda_t) + 1/(alpha_b_hx * (d_t/2)) ) * (T_hx[k] - T_w[k]);
    end if;
    if k == 1 then
      C_w[k]*der(T_w[k]) = N_hx*Qdot_hx[k] - Qdot_w_c[k] + Qdot_m[1];
    elseif k == N_w then
      C_w[k]*der(T_w[k]) = N_hx*Qdot_hx[k] - Qdot_w_c[k] - Qdot_m[N_w-1];
    else
      C_w[k]*der(T_w[k]) = N_hx*Qdot_hx[k] - Qdot_w_c[k] - Qdot_m[k - 1] + Qdot_m[k];
    end if;
    m_c[k]*c_c*der(T_c[k]) = (Qdot_w_c[k]) - Qdot_c_s[k];
// Im Kühlfall wird angenommen, dass der Massenstrom anderstherum in den Eisspeicher fließt (Testzwecke für Ergebnisse)
    if time < tN_h or time > t0_h then
      if k == 1 then
        m_hx/N_w*c_hx*der(T_hx[k]) = mdot_hx*c_hx*(inlet_hx.T - T_hx[k]) - Qdot_hx[k];
      else
        m_hx/N_w*c_hx*der(T_hx[k]) = mdot_hx*c_hx*(T_hx[k - 1] - T_hx[k]) - Qdot_hx[k];
      end if;
    else
      if k == N_w then
        m_hx/N_w*c_hx*der(T_hx[k]) = mdot_hx*c_hx*(inlet_hx.T - T_hx[k]) - Qdot_hx[k];
      else
        m_hx/N_w*c_hx*der(T_hx[k]) = mdot_hx*c_hx*(T_hx[k+1] - T_hx[k]) - Qdot_hx[k];
      end if;
    end if;
    V_s[k, 1]*c_s[k, 1]*der(T_s[k, 1]) = (Qdot_c_s[k]) - Qdot_s_s[k, 1];
    for l in 2:N_s - 1 loop
      V_s[k, l]*c_s[k, l]*der(T_s[k, l]) = (Qdot_s_s[k, l - 1]) - Qdot_s_s[k, l];
    end for;
    V_s[k, N_s]*c_s[k, N_s]*der(T_s[k, N_s]) = (Qdot_s_s[k, N_s - 1]) - Qdot_s_bd[k];
// Degree of icing
    p_ice[k] = max(0, -T_w[k]);
    der(x_ice[k]) = min(max(0, -T_w[k]),1) - x_ice[k];
  end for;
// Boundary condition
  T_s_bd = T_s_bd_min + ((-0.5*cos(2*Modelica.Constants.pi*((t0+time)/3600.0 - T_s_bd_shift)/8760.0)) + 0.5)*(T_s_bd_max - T_s_bd_min);
// Heat exchanger Outlet temperature
  if time < tN_h or time > t0_h then
    outlet_hx.T = T_hx[N_w];
  else
    outlet_hx.T = T_hx[1];
  end if;
// pressure
  delta_p = (8*l_t/N_hx*inlet_hx.mdot^2)/(rho_hx*Modelica.Constants.pi^2*d_t^5)*f_d;
  f_d = 0.3164 / ((Re+1e-9)^0.25);
  Re = (inlet_hx.mdot*d_t)/(mu*A_t);
  mu_lookup.u = inlet_hx.T;
  mu = mu_lookup.y[1];
  outlet_hx.p = inlet_hx.p;
// Monitoring
  Qdot_hx_sum = N_hx*sum(Qdot_hx);
  Qdot_s_s_sum = sum(Qdot_s_s);
  Qdot_c_s_sum = sum(Qdot_c_s);
  der(Q_hx) = N_hx*sum(Qdot_hx);
  der(Q_s) = sum(Qdot_c_s);
  der(Q_s_bd) = sum(Qdot_s_bd);
  Q_hx_kWh = to_kWh(Q_hx);
  Q_s_kWh = to_kWh(Q_s);
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-80, -100}, {80, 80}}), graphics = {Ellipse(lineThickness = 0.5, extent = {{-60, 68}, {60, 26}}), Line(points = {{-60, 48}, {-60, -80}, {60, -80}, {60, 46}}, thickness = 0.5), Line(points = {{20, -44}, {6, -38}, {52, -38}}, pattern = LinePattern.None, thickness = 0.5), Line(points = {{-86, 60}, {-76, 10}}, pattern = LinePattern.None, thickness = 0.5), Text(origin = {-4, -24}, extent = {{-36, 44}, {36, -44}}, textString = "E"), Line(origin = {2.96334, -61.6514}, rotation = 180, points = {{-36, 0}, {-14, 0}}, thickness = 0.5), Line(origin = {2.96334, -61.6514}, rotation = 180, points = {{-42, 0}, {-36, 6}, {-36, -6}, {-42, 0}}, thickness = 0.5), Line(origin = {9.64775, 17.8978}, rotation = 180, points = {{-36, 0}, {-14, 0}}, thickness = 0.5), Line(origin = {59.5171, 18.1794}, points = {{-42, 0}, {-36, 6}, {-36, -6}, {-42, 0}}, thickness = 0.5)}),
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-80, -100}, {80, 80}})));
end is;
