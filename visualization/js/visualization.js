//client.js : JavaScript code to that creates visualization.

// sets up line chart with scale, axes, and labels
function setup(margin, width, height) {

	// SVG
	var svg = d3.select("#plot").append("svg")
	    .attr("width", width + margin.left + margin.right)
	    .attr("height", height + margin.top + margin.bottom)
	  .append("g")
	    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	// Scales
	var x = d3.scale.linear()
	    .range([0, width]);
	var y = d3.scale.linear()
	    .range([height, 0]);
	var color = d3.scale.category10();

	// Axes
	var xAxis = d3.svg.axis()
	    .scale(x)
	    .orient("bottom")
	    .ticks(6);
	var yAxis = d3.svg.axis()
	    .scale(y)
	    .orient("left")
	    .ticks(6);

	// Graph
	svg.append("g")
	  .attr("class", "x axis")
	  .attr("transform", "translate(0," + height + ")")
	  .call(xAxis);
	svg.append("g")
	  .attr("class", "y axis")
	  .attr("transform", "translate(0,0)")
	  .call(yAxis);

	// Axes Labels
	svg.append("text")
	  .attr("class", "label")
	  .attr("text-anchor", "middle")
	  .attr("transform", "translate(" + width / 2 + "," + (height + margin.bottom / 2) + ")");              
	svg.append("text")
	  .attr("class", "label")
	  .attr("text-anchor", "middle")
	  .attr("transform", "translate(" + margin.left / -2 + "," + height / 2 + ")rotate(-90)");

	 // Tool Tip
	 var tip = d3.tip()
	 	.attr('class', 'd3-tip')
	 	.offset([-10, 0]);

	// return plot
	return {
		"svg": svg,
		"x": x,
		"y": y,
		"xAxis": xAxis,
		"yAxis": yAxis,
		"tip": tip
	};
}

// *************************************************************************************

var time = 300,
	color = ['#377eb8', '#e41a1c', '#4daf4a'];

// updates input plot labels and domain
function updatePlot (label, domain, plot) {
	// 7.1: Update Axes Labels
	plot.svg.selectAll("text.label")
	  .data(label)
	  .text(function(d) { return d; });

	// 7.2: Scale Domains
	plot.x.domain(domain.x);
	plot.y.domain(domain.y);

	// 7.3: Update Axis
	plot.svg.select('.x.axis')
	  .transition()
	  .duration(time)
	  .call(plot.xAxis);

	plot.svg.select(".y.axis")
	  .transition()
	  .duration(time)
	  .call(plot.yAxis);
}

// Adds data points to plot
function updatePoints (data, label, plot) {
		// update tool tip
		plot.tip.html(function(d,i) { return (label[0] + ":\t" + d.x.toString() + "<br>" + 
											  label[1] + ":\t" + d.y.toString()); });

	// JOIN new data with old elements.
	var alignment = plot.svg.selectAll("circle.stack")
	  .data(data[0]);

	// UPDATE old elements present in new data.
	alignment.transition(time)
	  .attr("cx", function(d){return plot.x(d.x)})
	  .attr("cy", function(d){return plot.y(d.y)});

	// ENTER new elements present in new data.
	alignment.enter().append("circle")
	  .attr("cx", function(d){return plot.x(d.x)})
	  .attr("cy", function(d){return plot.y(d.y)})
	  .attr("r", 3)
	  .attr("class", "stack")
	  .style("fill-opacity", 10^-6)
	  .on('mouseover', function(d){plot.tip.show(d,this)})
	  .on('mouseout', plot.tip.hide)
	.transition(time)
	  .style("fill-opacity", 1);

	// EXIT old elements not present in new data.
	alignment.exit()
	  .remove();

	plot.svg.call(plot.tip);
}

// updates core alignment paths
function updateAlignment (data, plot) {
	// DEFINE line function
	var lineFunction = d3.svg.line()
	  .x(function(d) { return plot.x(d.x); })
		  .y(function(d) { return plot.y(d.y); })
		  .interpolate("linear");

	// JOIN new data with old elements.
	var alignment = plot.svg.selectAll("path.stack")
	  .data(data);

	// UPDATE old elements present in new data.
	alignment.transition(time)
	  .attr("d", lineFunction);

	// ENTER new elements present in new data.
	alignment.enter().append("path")
	  .attr("d", lineFunction)
	  .attr("stroke", function(d,i) { return color[i]; })
	  .attr("stroke-width", 2)
	  .attr("stroke-opacity", 10^-6)
	  .attr("fill", "none")
	  .attr("class", "stack")
	.transition(time)
	  .style("stroke-opacity", 1);

	// EXIT old elements not present in new data.
	alignment.exit()
	  .remove();
}

// updates area stack
function updateStack (data, plot) {
	// DEFINE area function
		var areaFunction = d3.svg.area()
		  .x0(function(d) { return plot.x(d.x); })
	  .x1(function(d) { return plot.x(d.x); })
	  .y0(function(d) { return plot.y(d.y0); })
		  .y1(function(d) { return plot.y(d.y1); })
		  .interpolate("linear");

	// JOIN new data with old elements.
	var alignments = plot.svg.selectAll("path.stacks")
	  .data(data);

	// UPDATE old elements present in new data.
	alignments.transition(time)
	  .attr("d", areaFunction);

	// ENTER new elements present in new data.
	alignments.enter().append("path")
	  .attr("d", areaFunction)
	  .attr("fill-opacity", 10^-6)
	  .attr("fill", "black")
	  .attr("class", "stacks")
	.transition(time)
	  .style("fill-opacity", 0.25);

	// EXIT old elements not present in new data.
	alignments.exit()
	  .remove();
}

// updates confidence intervals
function updateInterval (data, plot) {
	// DEFINE variables
	var thick = 3;

	// JOIN new data with old elements.
	var intervals = plot.svg.selectAll("rect")
	  .data(data);

	// UPDATE old elements present in new data.
	intervals.transition(time)
	  .attr("x", function(d) { return plot.x(d.y0); })
	  .attr("y", function(d) { return plot.y(d.x) - thick/2; })
	  .attr("width", function(d) { return plot.x(d.y1) - plot.x(d.y0); });

	// ENTER new elements present in new data.
	intervals.enter().append("rect")
	  .attr("x", function(d) { return plot.x(d.y0); })
	  .attr("y", function(d) { return plot.y(d.x) - thick/2; })
	  .attr("width", function(d) { return plot.x(d.y1) - plot.x(d.y0); })
	  .attr("height", thick)
	  .attr("fill-opacity", 10^-6)
	  .attr("fill", "black")
	.transition(time)
	  .style("fill-opacity", 0.25);

	// EXIT old elements not present in new data.
	intervals.exit()
	  .remove();
}


// *************************************************************************************


// function that visualizes checked species
function visualize(stack) {
	// variables and data
	var curr = stack,
		ci = document.getElementById("confidence").checked,
		dp = document.getElementById("point").checked;

	// plot value plot
	var x_label = 'median',
		y_label = get_trait(),
		labels = [x_label,y_label],
		path = [extractPath(curr,x_label,y_label)],
		type = y_label == "delO18",
		area = ci && type ? [extractArea(ProbStack,"age","lower", "upper")] : [],
		bar = ci && !type ? extractArea(RadioCarbon[index],"depth","lower", "upper") : [],
		x_domain = [d3.min(path[0], function(d) { return d.x; }), 
			  		d3.max(path[0], function(d) { return d.x; })],
		y_domain = [d3.min(path[0], function(d) { return d.y; }), 
			  		d3.max(path[0], function(d) { return d.y; })];

	updatePlot(labels, {"x": x_domain, "y": y_domain}, value_plot);
	updateStack(area, value_plot);
	updateInterval(bar, value_plot);
	updateAlignment(path, value_plot);
	dp ? updatePoints(path, labels, value_plot) : updatePoints([[]], labels, value_plot);
	
	// plot timeline
	var labels = ["age","age offset"],
		path = [extractPath(curr,"median","lower","median"),
				extractPath(curr,"median","median","median"),
				extractPath(curr,"median","upper","median")],
		x_domain = [d3.min(path[0], function(d) { return d.x; }), 
			  		d3.max(path[0], function(d) { return d.x; })],
		y_domain = [-20,20];

	updatePlot(labels, {"x": x_domain, "y": y_domain}, time_plot);
	updateAlignment(path, time_plot);
}

// function returns array of x,y values
function extractPath(data,xVal,yVal0,yVal1){
	var points = [];
	for (var i = 0; i < data.length; i++) {
		var pt = data[i],
			x = pt[xVal],
			y = yVal1 ? pt[yVal0] - pt[yVal1] : pt[yVal0];
		if (!isNaN(x) && !isNaN(y)) {
			points.push({"x": x,"y": y});
		};
	};
	return points;
}

// function returns array of x,y values
function extractArea(data,xVal,yVal0,yVal1){
	var points = [];
	for (var i = 0; i < data.length; i++) {
		var pt = data[i],
			x = pt[xVal],
			y0 = pt[yVal0],
			y1 = pt[yVal1];
		if (!isNaN(x) && !isNaN(y0) && !isNaN(y1)) {
			points.push({"x": x,"y0": y0, "y1": y1});
		}
	};
	return points;
}

// function that returns random trait of iris dataset
function get_trait() {
	var e = document.getElementById("y-value");
	var value = e.options[e.selectedIndex].value;
	return value;
}

// 8b: Bind Visualize function to button on click
d3.select("#visualize").on("click", function() {
	var stack = getStack();
	visualize(stack);
});

// Gets Stack from selected radio buttons
function getStack() {
	var stacks = [SU81_18,MD03_2698,MD99_2334,MD95_2042,GeoB1032,GeoB1035,GeoB1214],
		index = document.getElementById("stack").selectedIndex;
	return run_hmm_stack(stacks[index], 0);
}


// 1: Create Plots
var margin = {top: 30, right: 30, bottom: 50, left: 60},
	width = 800 - margin.left - margin.right,
	height = 400 - margin.top - margin.bottom;
var value_plot = setup(margin, width, height);

var margin = {top: 10, right: 30, bottom: 60, left: 60},
	width = 800 - margin.left - margin.right,
	height = 200 - margin.top - margin.bottom;
var time_plot = setup(margin, width, height);


// TODO:
// Tool tip for value_plot
// D3 brushing for time plot
//	-https://bl.ocks.org/mbostock/34f08d5e11952a80609169b7917d4172
// General layout with bootstrap
// Write upload data function
// Contact Comp bio people about hosting
// Add resize methods
