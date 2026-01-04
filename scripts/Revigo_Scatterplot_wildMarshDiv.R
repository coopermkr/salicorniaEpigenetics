# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );

# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0002213","defense response to insect",0.15603406000624137,-6.048325595685048,3.5570418218417483,1.5563025007672873,-1.856985199745905,0.9335055241750008,0.35009449),
c("GO:0006518","peptide metabolic process",0.06687174000267487,-1.7374955785119484,-6.2863815366141536,1.2041199826559248,-2.2225731776106885,0.8903428360626815,0.09406659),
c("GO:0006575","modified amino acid metabolic process",0.5126833400205073,2.8767135762502294,-5.802824936847973,2.0644579892269186,-1.3458234581220394,0.8718951620501274,0.11322071),
c("GO:0006672","ceramide metabolic process",0.12036913200481478,-1.2045090551663777,-2.673372151213139,1.4471580313422192,-1.9665762445130504,0.7785568440244834,0.10972933),
c("GO:0006793","phosphorus metabolic process",3.3970843921358838,4.13243644913146,-3.0504160155825053,2.8825245379548803,-1.465973893943865,0.8483507949245882,0.13961335),
c("GO:0006797","polyphosphate metabolic process",0.004458116000178324,0.5961493460776055,7.821227181810362,0.3010299956639812,-3.3979400086720375,0.8464450000678954,0.32614192),
c("GO:0006799","polyphosphate biosynthetic process",0.004458116000178324,2.5397714319942493,6.390643371537692,0.3010299956639812,-3.3979400086720375,0.8092062986902813,0),
c("GO:0006817","phosphate ion transport",0.1827827560073113,-6.085315403230369,-0.9968314112789187,1.6232492903979006,-1.7878123955960423,0.9130000106246037,0.29885177),
c("GO:0008152","metabolic process",39.325041237573,-5.624656674291716,-4.011991431279596,3.9455670534423883,-2.167491087293764,0.9363495154756166,0.0261766),
c("GO:0009266","response to temperature stimulus",3.3168383041326734,-3.8895377937625946,6.225339284317604,2.8721562727482928,-1.47755576649368,0.9557956485912186,0.16121958),
c("GO:0009700","indole phytoalexin biosynthetic process",0.07132985600285319,4.548675061149958,3.931990489985377,1.2304489213782739,-2.1944991418415998,0.7248378358661876,0.16687265),
c("GO:0010193","response to ozone",0.15157594400606303,-4.586324959947948,3.786282966226812,1.5440680443502757,-1.869666231504994,0.9334882487677673,0.3524245),
c("GO:0016036","cellular response to phosphate starvation",0.4146047880165842,-5.394791290389724,4.988925931605694,1.9731278535996986,-1.4412914294668342,0.8943679119063243,0.22046999),
c("GO:0030643","intracellular phosphate ion homeostasis",0.022290580000891624,-3.891438897413744,-5.288783167909206,0.7781512503836436,-2.6989700043360187,0.977066106011776,-0),
c("GO:0032544","plastid translation",0.09362043600374481,3.9218189469917557,1.9010769731549009,1.3424226808222062,-2.0767559813697236,0.794174793642386,0.18633305),
c("GO:0032957","inositol trisphosphate metabolic process",0.031206812001248273,1.1484724670122384,7.644597719090767,0.9030899869919435,-2.5528419686577806,0.8092802922712874,0.36122349),
c("GO:0042398","modified amino acid biosynthetic process",0.12482724800499309,6.045207828512246,3.2375219722512574,1.462397997898956,-1.9507819773298183,0.7919134805158905,0.19396682),
c("GO:0042430","indole-containing compound metabolic process",0.37002362801480093,0.6894222105171323,-5.828904279577225,1.9242792860618816,-1.485452247339714,0.8752518872303011,0.10964638),
c("GO:0043043","peptide biosynthetic process",0.017832464000713297,5.420630882326528,5.102732820683368,0.6989700043360189,-2.795880017344075,0.8073193232567687,0.14206479),
c("GO:0043603","amide metabolic process",0.8827069680353082,5.233819725312987,-4.376430424154457,2.298853076409707,-2.5734887386354246,0.865889012262517,0.07679409),
c("GO:0044281","small molecule metabolic process",7.498551112299942,2.873025837707115,-2.886003380122389,3.226084115975824,-1.595166283380062,0.835962164121409,0.18962393),
c("GO:0045493","xylan catabolic process",0.09807855200392314,6.886964427311214,-1.17513455404199,1.3617278360175928,-2.0570004066339593,0.8433723255083727,0.09715713),
c("GO:0046686","response to cadmium ion",0.356649280014266,-3.4168143786162153,4.41328626213877,1.9084850188786497,-1.5003129173815961,0.9472751714413847,0.2687421),
c("GO:0052386","cell wall thickening",0.16495029200659803,0.6275849254945344,1.0218637063593143,1.5797835966168101,-1.8827287043442358,0.8821478476669068,0.26020134),
c("GO:0052542","defense response by callose deposition",0.12482724800499309,-5.542307699882564,2.1234973149821155,1.462397997898956,-1.9507819773298183,0.7819122991948748,-0),
c("GO:0072583","clathrin-dependent endocytosis",0.178324640007133,-6.744451002417073,-0.7262653533889056,1.6127838567197355,-1.7986028756795485,0.8417220249569142,0.2746206),
c("GO:1901659","glycosyl compound biosynthetic process",0.2541126120101645,5.622158206720742,2.2074456190142646,1.7634279935629373,-1.6536470255493614,0.7766217520853804,0.20926414));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$log_size <- as.numeric( as.character(one.data$log_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) );
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) ));
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);


# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("/path_to_your_file/revigo-plot.pdf");

