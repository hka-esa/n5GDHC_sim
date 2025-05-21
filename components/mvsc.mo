within n5GDHC.components;

model mvsc
  connectors.fluid_p inlet annotation(
    Placement(transformation(extent = {{-50, -10}, {-30, 10}}), iconTransformation(extent = {{-50, -10}, {-30, 10}})));
  connectors.fluid_p outlet_1 annotation(
    Placement(transformation(extent = {{-10, 30}, {10, 50}}), iconTransformation(extent = {{-10, 30}, {10, 50}})));
  connectors.fluid_p outlet_2 annotation(
    Placement(transformation(extent = {{-10, -50}, {10, -30}}), iconTransformation(extent = {{-10, -50}, {10, -30}})));
  Modelica.Blocks.Interfaces.RealInput pos annotation(
    Placement(visible = true, transformation(origin = {-52, 0}, extent = {{44, -18}, {84, 22}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{20, -10}, {0, 10}}, rotation = 0)));
equation
  assert(outlet_1.mdot >= (-1e-5), "MVC O1 negative massflow", level = AssertionLevel.error);
  assert(outlet_2.mdot >= (-1e-5), "MVC O2 negative massflow", level = AssertionLevel.error);
  assert(pos >= (-1e-5) and pos <= 1 + 1e-5, "MVC pos out of range", level = AssertionLevel.error);
  outlet_1.mdot = inlet.mdot*pos;
  outlet_1.T = inlet.T;
  outlet_1.p = inlet.p;
  outlet_2.mdot = inlet.mdot*(1 - pos);
  outlet_2.T = inlet.T;
  outlet_2.p = inlet.p;
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-60, -60}, {40, 60}}), graphics = {Line(points = {{-20, 40}, {20, -40}, {-20, -40}, {20, 40}, {-20, 40}}, thickness = 0.5), Line(points = {{0, 0}, {-40, 20}, {-40, -20}, {0, 0}}, thickness = 0.5), Line(points = {{-24, 6}, {-24, -6}, {-16, 0}, {-24, 6}}, thickness = 0.5), Text(extent = {{-16, 30}, {16, 12}}, textString = "1"), Text(extent = {{-16, -10}, {16, -28}}, textString = "0")}),
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-60, -60}, {40, 60}})));
end mvsc;
