// Initializes the PIXI application
function initialize_app(width, height, view, canvas=false){
	let app = new PIXI.Application(width, height, { 
		view: view,
		antialias: true, backgroundColor: 0xffffff, 
		resolution: devicePixelRatio, 
		autoResize: true, 
		forceCanvas: canvas,
		autoStart: true // note! 
	});
	return(app);	
}

// Create a viewport to handle panning, dragging, etc., add to the main stage
function initialize_vp(pn, app, width, height){
	var vp = pn.create_viewport(app, width, height);
	vp_params = {
		clampZoom: { minWidth: width/5, maxWidth: width*5, minHeight: height/5, maxHeight: height*5 }
	}
	// Zoom settings
	vp.drag().wheel(1e-12).clampZoom( vp_params.clampZoom ).decelerate();		
	return(vp);
}

HTMLWidgets.widget({
  name: 'pixiplex_vis',
  type: 'output',
  
  factory: function(pc_cont, width, height) {
		console.log(pc_cont);
		const pn = pixinet; console.log(pn);

    // Make the canvas view, append it to the widget container
  	var view = document.createElement('canvas');
	  view.width = width, view.height = height;
		view.onwheel = function(event){ event.preventDefault(); };
		pc_cont.appendChild(view); // append the canvas element

		// The rest of the global variables needed to be tracked
		var app = null, vp = null;
    var nodes = null, links = null, links_gfx = null; 
		var sim = null, sim_params = null; // the force simulation + parameters
		var ticker = null; 
		const in_shiny_mode = (HTMLWidgets.shinyMode && (typeof Shiny !== 'undefined'));
		// var Shiny = HTMLWidgets.shinyMode ? Shiny : null; 

		// Convert data.frame to array of objects
		const df_to_wide = (df) => {
			if (df.constructor === Array){ return df; } 
			else { return HTMLWidgets.dataframeToD3(df) }
		}

    // Return the factory API object 
    let initialized = false;
    return {
      enableForce: function() {
				console.log('enabling force');
				// Setup default force simulation
				if (!sim){
					sim = pn.enable_force();
					sim_params = pn.default_force_settings(); // Save default parameters as initialization
					pn.default_force_settings(sim, app, nodes, links); //  enact default settings 
				}

				// Register force additions to dragging callbacks for the nodes
				nodes.forEach((node) => { 
					pn.compose(pn.pixi_drag(node), pn.force_drag(sim))(pn.drag_dispatcher(node));
				});

				// Restart simulation
				ticker.on('tick', function(){
					links_gfx.clear();
					pn.draw_links(links, links_gfx);
					sim.tick(); // manually step through simulation
				})
				sim.alpha(1);
				ticker.call('animate');// start the animated rendering
				
				if (HTMLWidgets.shinyMode){
					console.log("in shiny mode");
				}
			},
			disableForce: function() {
				console.log('disabling force');
				// Setup default force simulation
				if (sim){
					ticker.stop();
					sim = null;
					// reassign dragging callbacks
					nodes.forEach((node) => { pn.compose(pn.pixi_drag(node))(pn.drag_dispatcher(node)); });
					// Restart animation loop w/o simulation
					ticker.on('tick', function(){
						links_gfx.clear();
						pn.draw_links(links, links_gfx);
					});
					ticker.call('animate');// start the animated rendering
				}
			},
			center: function(){
				vp.moveCenter(width/2, height/2);
			},
			setSimulation: function(x){
        if (sim){ pn.apply_sim(sim, x.params); }
			},
      setForce: function(x){
				if (sim && x.force){
					if (x.force.type == "forceLink"){ x.force.params.links = links; };
					pn.apply_force(sim, x.force);
				}
			},
			insertNodes: function(x){
				console.log(x)
				console.log(nodes)
				if ('nodes' in x){
					x.nodes = df_to_wide(x.nodes);					
					let ins_nodes = pn.differenceBy(x.nodes, nodes, 'id');
					ins_nodes = pn.generate_node_graphics(pn.scale_nodes(ins_nodes, width, height));
					ins_nodes.forEach((node) => {
						pn.compose(
							pn.pixi_drag(node), 
							sim ? pn.force_drag(sim) : pn.identity,
							in_shiny_mode ? pn.shiny_drag(Shiny) : pn.identity
						)(pn.drag_dispatcher(node));
					});
					pn.stage_items(vp, ins_nodes);
					nodes = pn.concat(nodes, ins_nodes);
					if (sim){ 
						sim.nodes(nodes);
						sim.alpha(1); 
					}
					console.log("here")
				}
			},
			// Removing nodes requires removing links attached to nodes as well
			removeNodes: function(x){
				if ('node_ids' in x){
					pn.remove_nodes(x.node_ids, nodes, links, vp);
					// Update simulation after any removal or addition
					if (sim){ 
						sim.nodes(nodes); 
						if (sim.force('link')){	sim.force('link').links(links); }
					}
				}
			},
			insertLinks: function(x){
				console.log(x.links)
				console.log(df_to_wide(x.links))
				if ('links' in x){
					x.links = df_to_wide(x.links);
					pn.resolve_links(nodes, x.links); // resolve source and target nodes before insertion 
					links = pn.unionWith(links, x.links, (o, p) => { return((o.source.id == p.source.id) && (o.target.id == p.target.id)) }); 
					if (sim && sim.force('link')){ 
						sim.force('link').links(links);
						sim.alpha(1);
					}
				}
			},
			removeLinks: function(x){
				console.log(x.links)
				console.log(df_to_wide(x.links))
				if ('links' in x){
					x.links = df_to_wide(x.links);
					pn.pullAllWith(links, x.links, (o, p) => { return((o.source.id == p.source) && (o.target.id == p.target)) });
					// TODO: check for all link forces!
					if (sim && sim.force('link')){ sim.force('link').links(links); }
					console.log("here")
				}
			},
			setNodeStyle: function(x){
				console.log(x);
				if (!('node_ids' in x)) { throw new Error("No node_ids given!") }
				if (!('style' in x)) { throw new Error("No node style supplied!") }
				if (x.node_ids.constructor == Number){ x.node_ids = [x.node_ids]; }
				let node_subset = pn.filter(nodes, (node) => { return pn.includes(x.node_ids, node.id) })
				x.style = df_to_wide(x.style);
				x.style = pn.map(x.style, (ns) => {
					let x_ns = {};
					if ("alpha" in ns){ x_ns.alpha = ns.alpha; }
					if ("color" in ns){ x_ns.color = ns.color; }
					if ("radius" in ns){ x_ns.radius = ns.radius; }
					if ("linesize" in ns){ 
						console.log("ls")
						if (!("lineStyle" in x_ns)){ x_ns.lineStyle = {}; };
						x_ns.lineStyle.size = ns.linesize; 
					}
					if ("linecolor" in ns){ 
						if (!("lineStyle" in x_ns)){ x_ns.lineStyle = {}; };
						x_ns.lineStyle.color = ns.linecolor; 
					}
					return x_ns;
				});
				let tmp = x.style
				console.log(tmp);
				console.log(node_subset);
				console.log(x.node_ids);
				pn.draw_nodes(node_subset, x.style);
			},
			setNodeText: function(x){
				if (!('node_ids' in x)) { throw new Error("No node_ids given!") }
				if (!('label' in x)) { throw new Error("No text labels supplied!") }
				if (!('text_style' in x)) { throw new Error("No text style supplied!") }
				let node_subset = pn.filter(nodes, (node) => { return pn.includes(x.node_ids, node.id) })
				x.text_style = df_to_wide(x.text_style);
				x.text_style = pn.map(x.text_style, (ts) => {
					let x_ts = {};
					if ("fontFamily" in ts){ x_ts.fontFamily = ts.fontFamily; }
					if ("fontSize" in ts){ x_ts.fontSize = ts.fontSize; }
					if ("fill" in ts){ x_ts.fill = ts.fill; }
					if ("align" in ts){ x_ts.align = ts.align;}
					return x_ts;
				});
				let tmp = x.text_style
				console.log(tmp);
				node_subset.forEach((node, i) => { 
					let text = new PIXI.Text(x.label[i], x.text_style);
					node.addChild(text);
				});
			},
			setLinkStyle: function(x){
				console.log(x);
				// if (!('source' in x) || !('target' in x)) { throw new Error("No source/targets given!") }
				if (!('style' in x)) { throw new Error("No link style supplied!") }
				let link_subset = pn.intersectionWith(links, x.links, (o, p) => { return((o.source.id == p.source) && (o.target.id == p.target)) })
				console.log(link_subset);
				pn.draw_links(link_subset, links_gfx, df_to_wide(x.style))
			},
			enableWeakTree: function(x){
				if (!('k' in x)) { throw new Error('Acceleration constant not given.'); }
				if (!('type' in x)) { throw new Error('Weak tree type not given'); }
				if (sim){
					let link_update = null;
					if (x.type == 0) { 
						link_update = function(link, k){ link.target.y += k; }
					} else if (x.type == 1){
						link_update = function(link, k){ link.source.y -= k; }
					} else {
						link_update = function(link, k){ link.source.y -= k; link.target.y += k; }
					}
					// Form weak tree by either pushing targets downward, sources upward, or both
					ticker.on('tick', function(e){
						links_gfx.clear();
						pn.draw_links(links, links_gfx);	
						links.forEach((link) => link_update(link, x.k * sim.alpha())); 
						sim.tick();
					});
				}
			},
      enableLasso: function(x){
        // Start lasso on double-click 
        let lasso_obj = pn.enable_lasso("pixir_container");
        lasso_obj.on("selected", function(nodes){
            // Highlight selected nodes 
            if (nodes) { 
              pn.draw_nodes(nodes, { color: 0xFFA500 });
              let node_container = pn.group_items(nodes);
              pn.stage_items(vp, [node_container]);
              let node_group = pn.enable_drag_container(node_container);
              console.log(node_group);
              pn.enable_drag(node_group).on("dragging", function(coords){
                console.log(this);
                this.x = coords.x; this.y = coords.y; 
                this.updateTransform();
              });
            };
        });
			},
      renderValue: function(x) {
				console.log("renderValue");
				console.log(x);

				// Convert to wide format
				if ('nodes' in x.plex) { x.plex.nodes = df_to_wide(x.plex.nodes) }
				if ('links' in x.plex) { x.plex.links = df_to_wide(x.plex.links) }
				const graph = x.plex;
				console.log(graph);

        // Initialization code
        if (!initialized){
					app = initialize_app(width, height, view, x.canvas);
					vp = initialize_vp(pn, app, width, height);
					app.stage.addChild(vp);// Enable dynamic resizing
					links_gfx = pn.generate_links_graphic();
					// Add drawing/clearing links to the ticker
					ticker = pn.ticker(app, vp); // the (pixi) simulation ticker
					ticker.on('tick', function(){
						links_gfx.clear();
						pn.draw_links(links, links_gfx);
					})
					document.getElementById(pc_cont.id).widget = this; // Attach widget to container to allow API calls
					initialized = true;  
				} // END INITIALIZATION 
				
				if (nodes || nodes) { pn.clear_stage(vp) }

				// Create interactive node + link graphics; scale to viewport
				links = graph.links;
				nodes = pn.generate_node_graphics(pn.scale_nodes(graph.nodes, width, height));
				pn.enable_interactive(nodes);
				
				// Resolve link source/targets + Enable node dragging
				pn.resolve_links(nodes, links);
				if (HTMLWidgets.shinyMode){ console.log("shiny mode"); }
				nodes.forEach((node) => {
					let dragger = (HTMLWidgets.shinyMode ? pn.compose(pn.pixi_drag(node), pn.shiny_drag(Shiny)) : pn.compose(pn.pixi_drag(node)));
					dragger(pn.drag_dispatcher(node));
				});

				// Add the graphics to the viewport/stage to be rendered + start the rendering
				pn.stage_items(vp, [links_gfx]); // add links to stage
				pn.stage_items(vp, nodes); // add nodes to stage
				ticker.call("animate");

				console.log("id: "+pc_cont.id);
				// Register Shiny handlers 
				if (HTMLWidgets.shinyMode) {
					var fxns = [
						'enableForce', 'disableForce', 'center', 'setSimulation', 'setForce',
						'insertNodes', 'removeNodes', 'insertLinks', 'removeLinks',
						'setNodeStyle', 'setLinkStyle', 'enableWeakTree', 'enableLasso', 
						'resize'
					];
					var addShinyHandler = function(fxn) {
						return function() {
							console.log("Registering shiny handler: ");
							console.log(pc_cont.id + ":"+fxn)
							Shiny.addCustomMessageHandler(
								pc_cont.id+":"+fxn, function(message) { // "pixiplex_vis:"
									var el = document.getElementById(message.id);
									if (el.widget) {
										console.log(message);
										el.widget[fxn](message);
									}
								}
							);
						}
					};
					for (var i = 0; i < fxns.length; i++) {
						addShinyHandler(fxns[i])();
					}
				} // end shiny handlers
				
				// Boilerplate: Call API methods passed in to enable chaining API calls
				if (x.hasOwnProperty('api')){
					var numApiCalls = x['api'].length;
					for (var i = 0; i < numApiCalls; i++) {
						var call = x['api'][i];
						var method = call.method;
						console.log("calling method: "+call.method);
						try { this[method](call); } catch(err) {}
					}
				} // api check
      },
      resize: function(width, height) { 
				if (vp){
					resize_f = pn.enable_resize(app, vp);
					scale_f = pn.make_scale(width, height);
				}
      }
		}; // return
	} // factory: function(pc_cont, width, height) 
}); // HTMLWidgets.widget

// window.pn = pn; 
// 				window.links = links; 
// 				window.nodes = nodes; 
// 				window.node_ids = x.node_ids;
// 				window.vp = vp;
// x.style = pn.map(x.style, (ns) => {
// 	x_ns = {
// 		color: ns.color, 
// 		lineStyle: { size: ns.linesize, color: ns.linecolor },
// 		radius: ns.radius,
// 		alpha: ns.alpha
// 	}
// 	x_ns.lineStyle = pn.clean(x_ns.lineStyle)
// 	return pn.clean(x_ns)
// });