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
revigo.data <- rbind(c("GO:0002239","response to oomycetes",0.45918594801836743,5.984335906007828,-2.2338687873392087,2.0170333392987803,-1.3555614105321614,0.9700904753279607,0.12674987),
c("GO:0006099","tricarboxylic acid cycle",0.24965449600998618,-2.3921495509811677,-7.268661866584262,1.7558748556724915,-1.6161846340195687,0.8639867391533557,0.1503925),
c("GO:0006520","amino acid metabolic process",1.894699300075788,-3.370774271629205,-5.260018762170931,2.629409599102719,-1.8477116556169435,0.8361706071811332,0.11244412),
c("GO:0008152","metabolic process",39.325041237573,7.220965486377663,0.6178840085615793,3.9455670534423883,-1.6143937264016879,0.9363052242550755,0.0398694),
c("GO:0009056","catabolic process",8.336676920333467,-3.6282616099366733,-3.2481953886638744,3.27207378750001,-1.3665315444204136,0.8307065442620705,0.22316787),
c("GO:0009308","amine metabolic process",0.494850876019794,-5.641439950935437,-5.175696021534799,2.0492180226701815,-1.324221658325915,0.873871540397252,0.12211281),
c("GO:0009309","amine biosynthetic process",0.21398956800855956,-2.1863276452732547,0.6225113034030324,1.6901960800285136,-1.6819366650372385,0.7796347490186476,0.09227854),
c("GO:0009987","cellular process",65.77950158263118,6.29607144448844,2.422257293049874,4.168968646554766,-1.3535962737769305,1,-0),
c("GO:0010600","regulation of auxin biosynthetic process",0.07578797200303151,5.604187409474337,4.390507313454771,1.255272505103306,-2.1307682802690238,0.8995346631640059,-0),
c("GO:0015698","inorganic anion transport",0.6553430520262137,-3.2773432538882212,5.81505541555572,2.1702617153949575,-6.436518914605589,0.9356139353151575,0),
c("GO:0015977","carbon fixation",0.12036913200481478,-6.501327746290358,2.7491779432481396,1.4471580313422192,-1.9318141382538383,0.8881362162451953,0.08725925),
c("GO:0022622","root system development",2.7506575721100264,2.245985535470578,5.746621550697216,2.790988475088816,-2.673664139071249,0.8832937440023706,0.49699056),
c("GO:0033512","L-lysine catabolic process to acetyl-CoA via saccharopine",0.013374348000534973,-7.151739514841789,-0.6067006309433335,0.6020599913279624,-2.8827287043442356,0.5596253630385002,0.38131679),
c("GO:0042128","nitrate assimilation",0.19169898800766796,-6.681927381189115,0.23282527689676596,1.6434526764861874,-3.826813731587726,0.6496172726319529,0.02090941),
c("GO:0042221","response to chemical",14.439837724577592,5.40824733895698,-3.3133077376797395,3.510545010206612,-1.3645162531850878,0.9572839866138075,0.22231597),
c("GO:0042430","indole-containing compound metabolic process",0.37002362801480093,-0.1555117638369115,-5.444013538027958,1.9242792860618816,-1.4473317838878068,0.8770969251885786,0.11840186),
c("GO:0042744","hydrogen peroxide catabolic process",0.4502697160180108,-6.560124408238252,-2.237123705595899,2.0086001717619175,-1.3645162531850878,0.8051603861757688,0.12088233),
c("GO:0044281","small molecule metabolic process",7.498551112299942,-2.680559010602091,-3.44562331786708,3.226084115975824,-2.3467874862246565,0.8328395468333979,0.17273807),
c("GO:0048527","lateral root development",0.5795550800231821,2.660238066920819,5.914361669103128,2.1172712956557644,-4.655607726314889,0.8458008638109127,-0),
c("GO:0051179","localization",12.861664660514466,3.093883894902587,1.7543856299271998,3.4602963267574753,-1.5287082889410615,1,-0),
c("GO:0055085","transmembrane transport",6.473184432258927,-3.066329489791748,5.3567521423089675,3.1622656142980214,-7.467245621007502,0.8716043599569872,0.48622625),
c("GO:0071941","nitrogen cycle metabolic process",0.22290580000891624,0.7422271537327791,-7.576317181841028,1.7075701760979363,-3.692503962086787,0.88233965075769,0.09258649),
c("GO:0072348","sulfur compound transport",0.19615710400784628,-3.7466979097109094,5.221477169567833,1.6532125137753437,-3.785156151952302,0.9424712726459502,0.33526483),
c("GO:1901698","response to nitrogen compound",0.8515001560340599,5.106618036378042,-4.573005009230223,2.2833012287035497,-2.5228787452803374,0.9453890711316447,0.27100724),
c("GO:1901701","cellular response to oxygen-containing compound",3.3658775801346352,4.705744468466599,-4.22408078588955,2.8785217955012063,-1.3861581781239307,0.8847530151333522,0.4352671),
c("GO:1902170","cellular response to reactive nitrogen species",0.057955508002318225,4.273957371242223,-5.257570941673684,1.146128035678238,-4.866461091629782,0.8916105011768342,-0),
c("GO:1902358","sulfate transmembrane transport",0.07578797200303151,-2.44459155436155,5.792765340500234,1.255272505103306,-4.6252516539898965,0.8967938419365963,0.3091586),
c("GO:1905393","plant organ formation",0.5171414560206856,1.9325139996052414,6.281428308539775,2.0681858617461617,-1.3053948010664311,0.8972149603680218,0.32898839),
c("GO:2001057","reactive nitrogen species metabolic process",0.21844768400873787,0.46677653320426643,-0.8962173729447216,1.6989700043360187,-3.709965388637482,0.8825393004420944,0.09145883));

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

