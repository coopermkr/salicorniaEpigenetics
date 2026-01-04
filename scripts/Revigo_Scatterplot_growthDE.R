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
revigo.data <- rbind(c("GO:0006099","tricarboxylic acid cycle",0.24965449600998618,5.02186235292229,3.5855767627977784,1.7558748556724915,-1.739928612014925,0.8885657170749304,0.10043649),
c("GO:0007623","circadian rhythm",0.5572645000222907,5.570780270042723,0.5915877195461106,2.100370545117563,-1.3990271043132518,1,-0),
c("GO:0009624","response to nematode",0.39231420801569256,0.4020431120457964,-6.086910739272,1.9493900066449128,-1.5451551399914898,0.9302068328267902,0.12273643),
c("GO:0009736","cytokinin-activated signaling pathway",0.38339797601533593,-2.5829696015237156,-6.6942949871579245,1.9395192526186185,-1.5559552040819238,0.8435966795404322,0.31572131),
c("GO:0010161","red light signaling pathway",0.04458116000178325,-1.931493402627612,-6.709653769301272,1.0413926851582251,-2.485452247339714,0.7867073028250421,0.01723293),
c("GO:0015695","organic cation transport",0.13820159600552806,-7.383762432798433,0.4209379888225354,1.505149978319906,-1.9956786262173574,0.865073780112883,0.32515367),
c("GO:0015698","inorganic anion transport",0.6553430520262137,-6.31442859294245,0.5846096343786538,2.1702617153949575,-3.0065637695023884,0.8468796314409758,0.27070924),
c("GO:0015837","amine transport",0.013374348000534973,-6.897622352508992,-1.2692999515102195,0.6020599913279624,-3.0078885122130505,0.8679238921875739,0.32532322),
c("GO:0015843","methylammonium transport",0.013374348000534973,-6.025399119029319,-1.5334397596568556,0.6020599913279624,-3.0078885122130505,0.8679238921875739,-0),
c("GO:0015977","carbon fixation",0.12036913200481478,3.039020915204916,5.972487112775127,1.4471580313422192,-2.0550240915879523,0.909073820600898,0.07240589),
c("GO:0019755","one-carbon compound transport",0.04012304400160492,-4.954372195016734,1.0580810065953106,1,-2.5316526695878427,0.8767487957777953,0.29385582),
c("GO:0033512","L-lysine catabolic process to acetyl-CoA via saccharopine",0.013374348000534973,-1.7457785362840283,7.290953641286618,0.6020599913279624,-3.0078885122130505,0.5559959330893272,0),
c("GO:0042128","nitrate assimilation",0.19169898800766796,-3.309471263598695,7.4850137179423175,1.6434526764861874,-1.8632794328435933,0.6853552764492542,0.38131679),
c("GO:0042744","hydrogen peroxide catabolic process",0.4502697160180108,-0.13786028771040093,7.012819617486085,2.0086001717619175,-1.486782399932061,0.7858200190247044,0.31801068),
c("GO:0048511","rhythmic process",0.6508849360260355,3.3702423754441773,-4.967757948918498,2.167317334748176,-1.332547047110046,1,-0),
c("GO:0071941","nitrogen cycle metabolic process",0.22290580000891624,5.2084430352982505,-2.389539560644264,1.7075701760979363,-1.7986028756795485,0.9046965123522359,0.09258649),
c("GO:0072348","sulfur compound transport",0.19615710400784628,-7.035087641779306,1.539836174422517,1.6532125137753437,-1.8446639625349381,0.8613589615136275,0.33526483),
c("GO:0072488","ammonium transmembrane transport",0.026748696001069945,-5.71988238553534,-0.5952149576297215,0.8450980400142568,-2.707743928643524,0.8353577039793246,0.33847477),
c("GO:0097272","ammonium homeostasis",0.026748696001069945,2.097923088482667,2.48907236306568,0.8450980400142568,-2.707743928643524,1,-0),
c("GO:1902358","sulfate transmembrane transport",0.07578797200303151,-6.105572133114551,1.912262282762977,1.255272505103306,-2.2549252084179425,0.8107218715497834,0.3091586),
c("GO:2001057","reactive nitrogen species metabolic process",0.21844768400873787,2.1395719154736583,-1.6920497833160062,1.6989700043360187,-1.8068754016455384,0.9048465938662578,0.08814633));

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

