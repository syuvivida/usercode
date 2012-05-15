import ROOT         # allows access to ROOT colors (e.g. ROOT.kRed)

##############################################################################
######## About Config Files ##################################################

## This file can be generated by running 'rootplot --config'

## Options are loaded in the following order:
##   1. from the command line
##   2. from the default configuration file
##   3. from ~/.rootplotrc
##   4. from configuration files specified on the command-line
## This leads to two major points worth understanding:
##   1. you may delete any lines you like in this file and they will still
##      be loaded correctly from the default
##   2. values specified here will superceed the same options from the
##      command-line
## Therefore, you could set, for example, 'xerr = True' in this file,
## and x-errorbars will always be drawn, regardless of whether '--xerr' is
## given on the command-line or not.  You can do this with any of the command-
## line options, but note that dashes are translated to underscores, so
## '--ratio-split=1' becomes 'ratio_split = 1'.

## Most global style options like default line widths can be set through
## a rootlogon.C, as described at:
##    http://root.cern.ch/drupal/content/how-create-or-modify-style

##############################################################################
######## Specifying Files and Targets ########################################

## You can specify the files to run on through the 'filenames' variable rather
## than entering them at the command-line, for example:
## filenames = ['histTTbar.root', 'histZmumu.root']

## Likewise, you can specify target histograms or directories here rather than 
## on the command-line, for example:
## targets = ['barrel/15to20', 'barrel/20to30']

## You might also want to specify fancy labels for the legend here rather 
## than on the command-line:
## legend_entries = [r'#bar{t}t', r'Z#rightarrow#mu#mu']

##############################################################################
######## Different Options for Different Targets #############################

## Leave these lists empty to have them automatically filled according to the
## command-line options.  Any list that is filled must be at least as long
## as the number of targets or it will throw an error.

line_colors = []                # normally filled by options.colors
fill_colors = []                # normally filled by options.colors
marker_colors = []              # normally filled by options.colors
marker_sizes = []         # in pixels
line_styles = []          # 1 (solid), 2 (dashed), 3 (dotted), 4 (dashdot), ...
fill_styles = []          # 0 (hollow), 1001 (solid), 2001 (hatched), ...
draw_commands = []        # a TH1::Draw option, include 'stack' to make stacked

##############################################################################
######## Global Style Options ################################################

## Colors can be specified as (r, g, b) tuples (with range 0. to 1. or range
## 0 to 255), or ROOT color constants (ROOT.kBlue or 600)

colors = [
    ## a default set of contrasting colors the author happens to like
#    (0,0,0), #kBlack for data
#    ROOT.kViolet-5, # kPink+3 for qcd 
#    ROOT.kOrange+7, # kYellow-2 for tt
#    ROOT.kOrange-2, # kTeal-5 for zjets
#    ( 0,51,153), # kAzure-6
#    (145,  83, 207), # purple
#    (231, 139,  77), # orange
#    (114, 173, 117), # green
#    ( 67,  77,  83), # dark grey
    ROOT.kBlack,
#    ( 67,  77,  83), # dark grey
    ROOT.kRed,
    ROOT.kGreen,
    ROOT.kBlue,
    ROOT.kOrange,
    (145,  83, 207),
    ROOT.kCyan,
    ( 67,  77,  80), # dark grey
    ]

## Used when --marker_styles is specified; more info available at:
## http://root.cern.ch/root/html/TAttMarker.html
marker_styles = [
     20, # circle
     25, # circle
     2, # circle
    26, # triangle
     5, # x
    30, # five-pointed star
    27, # diamond
    28, # cross
     3, # asterisk
    ]

#### Styles for --data
data_color = (0,0,0)      # black
mc_color = (ROOT.kOrange-2) # used when there are exactly 2 targets; set to
                          # None to pick up the normal color
data_marker = 23           # marker style (circle)

#### Settings for --ratio-split or --efficiency-split
ratio_max = 1.5
ratio_min = 0.5
ratio_fraction = 0.50  # Fraction of the canvas that bottom plot occupies
#ratio_label = '    Data/MC' # Label for the bottom plot
ratio_label = 'Ratio to %(ratio_file)s' # Label for the bottom plot
efficiency_label = 'Efficiency vs. %(ratio_file)s'

#### Titles produced by --area-normalize and --normalize
area_normalized_title = 'Fraction of Events in Bin'
target_normalized_title = 'Events Normalized to %(norm_file)s'

#### Overflow and underflow text labels
overflow_text = ' Overflow'
underflow_text = ' Underflow'

#### Define how much headroom to add to the plot
top_padding_factor = 1.2
top_padding_factor_log = 5.    # used when --logy is set

#### Plotting options based on histogram names
## Apply options to histograms whose names match regular expressions
## The tuples are of the form (option_name, value_to_apply, list_of_regexs)
## ex: to rebin by 4 all histograms containing 'pt' or starting with 'eta':
##    ('rebin', 4, ['.*pt.*', 'eta.*'])
options_by_histname = [
    ('area_normalize', True, []),
                       ]

#### Legend
legend_width = 0.26 #38        # Fraction of canvas width
legend_entry_height = 0.03 # Fraction of canvas height
max_legend_height = 0.4    # Fraction of canvas height
legend_left_bound = 0.20   # For left justification
legend_right_bound = 0.94  # For right justification
legend_upper_bound = 0.93  # For top justification
legend_lower_bound = 0.15  # For bottom justification
legend_codes = { 1 : 'upper right',
                 2 : 'upper left',
                 3 : 'lower left',
                 4 : 'lower right',
                 5 : 'right',
                 6 : 'center left',
                 7 : 'center right',
                 8 : 'lower center',
                 9 : 'upper center',
                10 : 'center',
                }

#### Page numbers
numbering_size_root = 0.03  # Fraction of canvas width
numbering_align_root = 33   # Right-top adjusted
numbering_x_root = 0.97     # Fraction of canvas width
numbering_y_root = 0.985    # Fraction of canvas height

#### Draw style for TGraph
draw_graph = 'ap'

#### This code snippet will be executed after the histograms have all
#### been drawn, allowing you to add decorations to the canvas
decoration_root = '''
## Draw a line to indicate a cut
#line = ROOT.TLine(5.,0.,5.,9.e9)
#line.Draw()
## Add a caption
##tt = ROOT.TText()
##tt.DrawTextNDC(0.6, 0.90, "2.050 /fb, 7 TeV")
latex = ROOT.TLatex()
##tt.SetNDC()
##tt.DrawLatex(0.6, 0.90, "2.050 /fb, #sqrt{s} = 7 TeV")
latex.SetNDC()
latex.SetTextSize(0.04)
latex.SetTextAlign(11)
##latex.DrawLatex(0.48,0.88,"#int #font[12]{L} dt = 2.050 fb^{-1}")
latex.SetTextAlign(31)
#latex.DrawLatex(0.94,0.96,"4.890 fb^{-1}  at   #sqrt{s} = 7 TeV")
latex.DrawLatex(0.94,0.96,"MCFM with CT10as scale m(34)")
#latex.SetTextAlign(11)
#latex.DrawLatex(0.16,0.96,"CMS preliminary 2011")

'''

##############################################################################
######## HTML Output #########################################################

#### Number of columns for images in HTML output
ncolumns_html = 2

#### Provide a template for the html index files
html_template=r'''
    <html>
    <head>
    </head>
    <body>
      %(back_nav)s
      <ul>
          %(forward_nav)s
      </ul>
    <h1>Images</h1>
    %(plots)s
    </body>
    </html>
'''
