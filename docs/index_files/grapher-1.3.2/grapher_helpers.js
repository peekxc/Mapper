/**
 * Created by mpiekenbrock on 3/3/18.
 */

// Setup D3's force layout
function D3_forceNet(network){
	var force = d3.layout.force()
		.nodes(network.nodes)
		.links(network.links)
		.size([width, height])
		.on('tick', onTick)
		.charge(-5000)
		.gravity(0.005)
		.linkStrength(0.2)
		.linkDistance(100)
		.friction(0.02)
		.start();
	return(force);
}
