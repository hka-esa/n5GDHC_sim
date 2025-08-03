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

within n5GDHC.constructs;

/* 
The parameter data is taken from literature, public sources or data sheets in order to reflect the specific conditions at the site.
Sources:
[1] Kuchling, H.: Taschenbuch der Physik, Hanser München, 22, 2022.
[2] Dehner, U.: Boden und Energiewende, Springer Fachmedien Wiesbaden, 2015.
[3] Deutscher Wetterdienst: Open Data, https://opendata.dwd.de/, 2024.
*/

model hn_con_h
  parameter Integer N_s = 5 "Number of sorrounding soil layers";
  parameter Modelica.Units.SI.Length l_n = 20.0 "Length of network tubes";
  parameter Modelica.Units.SI.Length l_h = 5.0 "Length of household tube";
  parameter Modelica.Units.SI.Radius r_n = 0.06 "Radius of network tube";
  parameter Modelica.Units.SI.Radius r_h = 0.01 "Radius of household tube";
  parameter Modelica.Units.SI.Radius t_t_n = 0.005 "Thickness of network tube";
  parameter Modelica.Units.SI.Radius t_t_h = 0.005 "Thickness of household tube";
  parameter Modelica.Units.NonSI.Temperature_degC T_init = 10.0;
  parameter Modelica.Units.NonSI.Temperature_degC T_s_0 = 10.0;
  parameter Modelica.Units.NonSI.Temperature_degC T_s_bd_min = 5.0 "Minimum yearly soil temperature of nearest DWD weather station [3]";
  parameter Modelica.Units.NonSI.Temperature_degC T_s_bd_max = 20.0 "Maximum yearly soil temperature of nearest DWD weather station [3]";
  parameter Modelica.Units.SI.ThermalConductivity lambda_s = 1.2 "Therm. conductivity soil [1,2]";
  parameter Modelica.Units.SI.Time t0 = 0.0;
  parameter Real ratio_s_w = 0.10;
  components.hn_split hn_split annotation(
    Placement(visible = true, transformation(origin = {8, -14}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  components.hn_tubes hn_tubes(N_s = N_s, l = l_n, r_t = r_n, t_t = t_t_n, T_f_init = T_init, lambda_s = lambda_s, t0 = t0, ratio_s_w = ratio_s_w, T_s_0 = T_s_0, T_s_bd_min = T_s_bd_min, T_s_bd_max = T_s_bd_max) annotation(
    Placement(visible = true, transformation(origin = {-47, -31}, extent = {{-37, 35}, {37, -35}}, rotation = 0)));
  components.hn_tubes h_tubes(N_s = N_s, l = l_h, r_t = r_h, t_t = t_t_h, T_f_init = T_init, lambda_s = lambda_s, t0 = t0, ratio_s_w = ratio_s_w, T_s_0 = T_s_0, T_s_bd_min = T_s_bd_min, T_s_bd_max = T_s_bd_max) annotation(
    Placement(visible = true, transformation(origin = {24, 37}, extent = {{33, -32}, {-33, 32}}, rotation = -90)));
  components.hn_merge hn_merge annotation(
    Placement(visible = true, transformation(origin = {37, -48}, extent = {{-11, -12}, {11, 12}}, rotation = 0)));
  connectors.fluid_p inlet_hn_sup annotation(
    Placement(visible = true, transformation(origin = {-91, -13}, extent = {{-7, -7}, {7, 7}}, rotation = 0), iconTransformation(origin = {-60, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  connectors.fluid_p outlet_hn_ret annotation(
    Placement(visible = true, transformation(origin = {-93, -49}, extent = {{-7, -7}, {7, 7}}, rotation = 0), iconTransformation(origin = {-58, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  connectors.fluid_p outlet_h_sup annotation(
    Placement(visible = true, transformation(origin = {7, 79}, extent = {{-7, -7}, {7, 7}}, rotation = 0), iconTransformation(origin = {0, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  connectors.fluid_p inlet_h_ret annotation(
    Placement(visible = true, transformation(origin = {39, 79}, extent = {{-7, -7}, {7, 7}}, rotation = 0), iconTransformation(origin = {0, 82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  connectors.fluid_p outlet_hn_sup annotation(
    Placement(visible = true, transformation(origin = {65, -13}, extent = {{-7, -7}, {7, 7}}, rotation = 0), iconTransformation(origin = {60, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  connectors.fluid_p inlet_hn_ret annotation(
    Placement(visible = true, transformation(origin = {67, -49}, extent = {{-7, -7}, {7, 7}}, rotation = 0), iconTransformation(origin = {62, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(hn_tubes.outlet_hn_supply, hn_split.inlet) annotation(
    Line(points = {{-22, -14}, {0, -14}}, color = {28, 108, 200}));
  connect(hn_split.outlet_1, h_tubes.inlet_hn_supply) annotation(
    Line(points = {{8, -6}, {8, 15}}, color = {28, 108, 200}));
  connect(h_tubes.outlet_hn_supply, outlet_h_sup) annotation(
    Line(points = {{8, 59}, {8, 80}}, color = {28, 108, 200}));
  connect(inlet_h_ret, h_tubes.inlet_hn_return) annotation(
    Line(points = {{40, 80}, {40, 59}}));
  connect(h_tubes.outlet_hn_return, hn_merge.inlet_1) annotation(
    Line(points = {{40, 16}, {37, 16}, {37, -40}}, color = {28, 108, 200}));
  connect(hn_tubes.inlet_hn_supply, inlet_hn_sup) annotation(
    Line(points = {{-72, -14}, {-90, -14}, {-90, -12}}, color = {28, 108, 200}));
  connect(hn_split.outlet_2, outlet_hn_sup) annotation(
    Line(points = {{16, -14}, {66, -14}, {66, -12}}, color = {28, 108, 200}));
  connect(hn_merge.outlet, hn_tubes.inlet_hn_return) annotation(
    Line(points = {{30, -48}, {-22, -48}}, color = {28, 108, 200}));
  connect(inlet_hn_ret, hn_merge.inlet_2) annotation(
    Line(points = {{68, -48}, {44, -48}}));
  connect(hn_tubes.outlet_hn_return, outlet_hn_ret) annotation(
    Line(points = {{-72, -48}, {-92, -48}}, color = {28, 108, 200}));
  annotation(
    Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})),
    Icon(graphics = {Line(origin = {5.94641, 50.3156}, points = {{-42, 0}, {-36, 6}, {-36, -6}, {-42, 0}}, thickness = 0.5), Rectangle(origin = {10, 50}, lineThickness = 0.5, extent = {{-68, 30}, {52, -30}}), Line(origin = {42.0287, 50.3441}, points = {{-42, 0}, {-36, 6}, {-36, -6}, {-42, 0}}, thickness = 0.5), Line(origin = {0.0245658, 20.3246}, rotation = -90, points = {{-30, 0}, {-36, 6}, {-36, -6}, {-30, 0}}, thickness = 0.5), Line(origin = {0.0240458, 28.3627}, rotation = -90, points = {{-40, 0}, {-28, 0}}, thickness = 0.5), Line(origin = {-5.95219, -50.4347}, rotation = 180, points = {{-42, 0}, {-36, 6}, {-36, -6}, {-42, 0}}, thickness = 0.5), Line(origin = {-41.9683, -50.6038}, rotation = 180, points = {{-42, 0}, {-36, 6}, {-36, -6}, {-42, 0}}, thickness = 0.5), Line(origin = {0.0288379, -2.5147}, rotation = 90, points = {{-30, 0}, {-36, 6}, {-36, -6}, {-30, 0}}, thickness = 0.5), Line(origin = {-55.967, -50.6103}, rotation = 180, points = {{-50, 0}, {-20, 0}}, thickness = 0.5), Rectangle(origin = {-8, -50}, rotation = 180, lineThickness = 0.5, extent = {{-68, 30}, {52, -30}}), Line(origin = {0.0329359, -10.5696}, rotation = 90, points = {{-40, 0}, {-28, 0}}, thickness = 0.5), Ellipse(origin = {-8, -50}, rotation = 180, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{34, 30}, {68, -30}}), Ellipse(origin = {10, 50}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{34, 30}, {68, -30}}), Ellipse(origin = {-7, 50}, rotation = 180, fillColor = {255, 255, 255}, pattern = LinePattern.Dash, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{34, 30}, {68, -30}}, startAngle = 90, endAngle = 270, closure = EllipseClosure.Chord), Ellipse(origin = {-109, 50}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{34, 30}, {68, -30}}, startAngle = 90, endAngle = 270, closure = EllipseClosure.Chord), Ellipse(origin = {-20, 50}, rotation = 180, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{25, 29.75}, {49, -29.75}}), Ellipse(origin = {111, -50}, rotation = 180, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{34, 30}, {68, -30}}, startAngle = 90, endAngle = 270, closure = EllipseClosure.Chord), Ellipse(origin = {9, -50}, fillColor = {255, 255, 255}, pattern = LinePattern.Dash, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{34, 30}, {68, -30}}, startAngle = 90, endAngle = 270, closure = EllipseClosure.Chord), Ellipse(origin = {96, -50}, rotation = 180, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{25, 29.75}, {49, -29.75}}), Line(origin = {-19.9614, -50.6085}, rotation = 180, points = {{-50, 0}, {-20, 0}}, thickness = 0.5), Line(origin = {-49.9869, 50.3368}, rotation = 180, points = {{-50, 0}, {-20, 0}}, thickness = 0.5), Line(origin = {-13.9146, 50.3041}, rotation = 180, points = {{-50, 0}, {-20, 0}}, thickness = 0.5)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
end hn_con_h;
