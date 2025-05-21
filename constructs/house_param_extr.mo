within n5GDHC.constructs;

model house_param_extr
  parameter Real idx_r = 0;
  parameter Integer idx = integer(idx_r);
  parameter String path_Qdot = Modelica.Utilities.Files.loadResource("modelica://n5GDHC/input/Qdot_" + String(idx) + ".txt");
  parameter Modelica.Units.NonSI.Temperature_degC T_init = 10.0 "Initial fluid temperature";
  parameter Modelica.Units.NonSI.Temperature_degC dT = 1.0 "Temperature difference for hx";
  components.pump pump_hn annotation(
    Placement(visible = true, transformation(origin = {-7, -29}, extent = {{9, -9}, {-9, 9}}, rotation = 90)));
  connectors.fluid_p inlet_hn annotation(
    Placement(visible = false, transformation(origin = {-26, 44}, extent = {{16, -96}, {22, -90}}, rotation = 0), iconTransformation(origin = {0, -2}, extent = {{20, -68}, {40, -48}}, rotation = 0)));
  connectors.fluid_p outlet_hn annotation(
    Placement(visible = false, transformation(origin = {-30, 44}, extent = {{38, -96}, {44, -90}}, rotation = 0), iconTransformation(origin = {0, -2}, extent = {{-40, -68}, {-20, -48}}, rotation = 0)));
  components.hx_param_qdot hx_param_qdot(T_init = T_init) annotation(
    Placement(visible = true, transformation(origin = {2, 3}, extent = {{-24, -13}, {24, 13}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable Qdot_ext(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic,fileName = path_Qdot, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, tableName = "Qdot", tableOnFile = true) annotation(
    Placement(visible = true, transformation(origin = {-30, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  inlet_hn.mdot = outlet_hn.mdot;
  pump_hn.mdot = abs(Qdot_ext.y[1]) / (hx_param_qdot.c_B * dT);
  connect(pump_hn.inlet, inlet_hn) annotation(
    Line(points = {{-7, -36}, {-7, -49}}, color = {28, 108, 200}));
  connect(Qdot_ext.y[1], hx_param_qdot.Qdot_ext) annotation(
    Line(points = {{-19, 30}, {1, 30}, {1, 11}}, color = {0, 0, 127}));
  connect(pump_hn.outlet, hx_param_qdot.inlet_B) annotation(
    Line(points = {{-7, -22}, {-7, -6}}, color = {28, 108, 200}));
  connect(hx_param_qdot.outlet_B, outlet_hn) annotation(
    Line(points = {{11, -6}, {11, -49}}, color = {28, 108, 200}));
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-80, -80}, {80, 80}}), graphics = {Rectangle(lineThickness = 0.5, extent = {{-60, 40}, {60, -60}}), Line(points = {{-60, 40}, {0, 60}, {60, 40}}, thickness = 0.5), Rectangle(origin = {0, -35}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-44, 19}, {44, -19}}), Line(origin = {16.6185, -55.8114}, points = {{-60, 40}, {28, 2}}, thickness = 0.5), Text(origin = {1, -45}, extent = {{-40, 7}, {-19, -5}}, textString = "B"), Text(origin = {61, -27}, extent = {{-40, 7}, {-19, -5}}, textString = "A")}),
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-80, -80}, {80, 80}}), graphics = {Text(origin = {60, 52},extent = {{-88, -44}, {-80, -52}}, textString = "HX"), Text(origin = {-26, 34}, extent = {{44, -60}, {78, -92}}, textString = "Heating
network")}));
end house_param_extr;