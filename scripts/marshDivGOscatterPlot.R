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
revigo.data <- rbind(c("GO:0002213","defense response to insect",0.0005138624181922514,6.867760389800543,-1.950890035348646,2.2600713879850747,-1.8181564120552274,0.9038058876869393,0.33895024),
c("GO:0006518","peptide metabolic process",0.38701789342754633,6.52269395068341,4.834377721464634,5.134565949067227,-2.185086818724926,0.9527783529335364,0.05818439),
c("GO:0006575","modified amino acid metabolic process",1.0391348532834788,4.309849042138549,-7.045987409797914,5.563504816684922,-2.9507819773298185,0.9485021953906324,0.06490785),
c("GO:0006643","membrane lipid metabolic process",0.8791731732247813,6.267016361686911,-5.443878086584604,5.490907548684358,-1.3053948010664311,0.9230564093853658,0.07207014),
c("GO:0006749","glutathione metabolic process",0.3259136899648291,-0.46297270124876794,-0.5304846034246018,5.0599381049925,-3.486782399932061,0.8447354879466887,0.05715436),
c("GO:0006790","sulfur compound metabolic process",2.195890258974167,0.25133847707195495,7.439503789160145,5.888442912219316,-1.9665762445130504,0.9447099425604463,0.07113475),
c("GO:0006793","phosphorus metabolic process",8.027757428321303,-6.136625631193294,-2.2817755170469183,6.4514261132466615,-1.3946949538588904,0.936696498207186,0.12455894),
c("GO:0006797","polyphosphate metabolic process",0.05285969262056093,-2.437330796637594,5.621691242689403,4.269979676645324,-3.359518563029578,0.9191044363845651,0.04813428),
c("GO:0006799","polyphosphate biosynthetic process",0.0298977078783569,-3.4658708760501256,4.988077388091639,4.02251085043403,-3.359518563029578,0.8972248205067511,0.29626057),
c("GO:0006817","phosphate ion transport",0.23164690690595885,1.2663980358894678,-6.557324114584249,4.911663546757045,-3.8477116556169433,0.9522919014336517,-0),
c("GO:0006898","receptor-mediated endocytosis",0.20247030982440214,0.7950627114973948,-5.444005847959121,4.853199155575347,-1.7304870557820837,0.920421961298327,0.21799363),
c("GO:0008152","metabolic process",54.415440901336254,-5.2947101454216945,-4.35589839180727,7.282553893608798,-1.6143937264016879,0.9929076516861575,0.00771111),
c("GO:0009266","response to temperature stimulus",0.37155375756073317,6.891349022771279,1.306392589678986,5.1168566946868985,-1.4056074496245734,0.90966473108778,0.25402444),
c("GO:0009404","toxin metabolic process",0.11547141511123447,2.5982838122234,7.176572745422406,4.6093168842930865,-3.5016894462103996,0.8787497531770221,0.05163483),
c("GO:0009635","response to herbicide",0.005558798976908442,6.089139195296818,-0.5634261123357732,3.2920344359947364,-2.1567672219019904,0.8873137124073561,0.27678133),
c("GO:0009987","cellular process",87.14322191622645,4.800525962393403,6.2837482590415075,7.487065345620153,-1.3535962737769305,1,-0),
c("GO:0015698","inorganic anion transport",0.848466344964188,1.7978549107438588,-5.768676492351584,5.475467792116027,-2.7471469690201067,0.9481440033053921,0.24681557),
c("GO:0016036","cellular response to phosphate starvation",0.06901654909532391,7.068710277622748,-0.7793986986996354,4.385802823325019,-3.1450869776921446,0.8424901984366238,0.00314292),
c("GO:0019184","nonribosomal peptide biosynthetic process",0.1888458581950949,-6.02053277386224,1.9634682258494136,4.822945711717861,-2.8827287043442356,0.8843042205733929,0.05411224),
c("GO:0019748","secondary metabolic process",0.6785567427323104,-3.3637141489860958,-0.985719521369392,5.3784197059814165,-2.2932822176632413,0.9504392247019186,0.06182518),
c("GO:0031667","response to nutrient levels",0.528789979489792,7.429043562804884,0.6502177512232044,5.27011726695538,-2.1739251972991736,0.9072130824099812,0.21987267),
c("GO:0032957","inositol trisphosphate metabolic process",0.009337533112896767,-2.7800564869708007,6.141560362717752,3.5171958979499744,-2.515700160653214,0.9108898479031542,0.27576801),
c("GO:0033037","polysaccharide localization",0.0782632336038994,2.280074735084865,-4.44649064452624,4.440405260103954,-1.6736641390712486,0.9383782999892625,0.19214761),
c("GO:0042430","indole-containing compound metabolic process",0.32694993185786875,-1.8716974654061067,-7.2070512186272335,5.061316740851394,-1.4473317838878068,0.9534381213263345,0.05717306),
c("GO:0042435","indole-containing compound biosynthetic process",0.2244159258057618,-5.543805920181961,3.405766354657811,4.897890886286548,-1.6020599913279623,0.9096954418936364,0.14691498),
c("GO:0042592","homeostatic process",2.0952811077261173,-3.7060905206330093,-6.0084536705625995,5.868074604162045,-1.7569619513137056,1,-0),
c("GO:0043603","amide metabolic process",1.8676684466558002,-2.7643245253891218,-3.7634494113504644,5.818132160374245,-4.1249387366083,0.9455765754900828,-0),
c("GO:0044281","small molecule metabolic process",14.48979594154301,-6.558590247413091,-0.750464711535265,6.7078940730325085,-1.489454989793388,0.9323658910928989,0.09647906),
c("GO:0045493","xylan catabolic process",0.0931652437314747,3.55391438041644,4.406820983048681,4.516098877052355,-2.0190880622231564,0.9156344924016178,0.05062349),
c("GO:0046686","response to cadmium ion",0.03216324494861887,5.556702884013442,0.3866262478205818,4.054229909863397,-1.4634415574284698,0.9018708604976646,0.3189107),
c("GO:0051179","localization",20.31279969393105,0.34252292903892545,2.1390728862426864,6.854601565181338,-1.5287082889410615,1,-0),
c("GO:0052544","defense response by callose deposition in cell wall",0.0007523400045356166,4.940275921322165,-2.1430888775900803,2.424881636631067,-2.1062382379420566,0.7688223785179411,0.29975261),
c("GO:0055062","phosphate ion homeostasis",0.03825861849479988,0.9865391787589963,5.042954684666853,4.129593228367933,-4.480172006224281,0.8961908776809364,0),
c("GO:1901659","glycosyl compound biosynthetic process",0.19687176458310315,-4.371904346085709,2.2031704377636334,4.841021415257488,-1.6161846340195687,0.9028257427390805,0.14488711));

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

