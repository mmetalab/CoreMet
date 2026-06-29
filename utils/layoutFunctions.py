from cProfile import label
from pydoc import classname
from dash import dcc, html
import dash_bootstrap_components as dbc
from dash import Dash, dash_table, Input, Output, callback
import dash_bio as dashbio

# ------------------------------------------------------------------------------
# make_ FUNCTIONS
# These functions are callse only once to create the backbone of the graphical 
# plots and elements. Then each element is updated based on callbacks that call 
# "update_" functions. 
# ------------------------------------------------------------------------------

def make_Subtitle(string):
    """
    Makes a subtitle with a line underneath to divide sections
    """
    subtitle = html.Div([
        html.H2(string,
            className="mt-2 mb-0 fw-bolder"),
        html.Hr(className="my-0")   
    ])
    return subtitle

def make_NavBar():
    """
    Makes the navigation bar
    """
    navbar = dbc.NavbarSimple(
        children=[
            dbc.NavItem(dbc.NavLink('Introduction', href='/introduction')),
            dbc.NavItem(dbc.NavLink('Database', href='/database')),
            dbc.NavItem(dbc.NavLink('Interactions', href='/interactions')),
            dbc.NavItem(dbc.NavLink('Docking', href='/docking')),
            dbc.DropdownMenu(
                children=[
                    dbc.DropdownMenuItem('Cite', id='citeDropdown'),
                    dbc.DropdownMenuItem('About Us', id='aboutUsDropdown'),
                ],
                nav=True,
                in_navbar=True,
                label='More',
            ),
        ],
        brand='CoreMet',
        brand_href='/introduction',
        color='primary',
        fixed='top',
        dark=True,
        style={'height':'40px'},
        className=''
    )
    return navbar

def make_IntroHeader(idFunc):
    """
    Makes the header for the WFA page
    """
    # Main Title
    header = html.Div(
        dbc.Container([
            html.Div([
                html.H2("Enzymatic reaction link learning", className="display-4"),
                dbc.Button("Cite", id=idFunc("btn_citeHeader"),color="success", outline=True)
            ], className='d-flex justify-content-between align-items-center mb-0'),
            
            html.Hr(className="mt-0 mb-1"),
            html.Div([
                html.P("Enzymatic reaction link prediction by varaitional graph autoencoders"),
                html.H4(id=idFunc('moreInfoIcon'), className="fa-solid fa-book-open ms-3 mt-1 primary")
            ], className='d-flex mb-0'),
            dbc.Collapse(
                html.Div([
                    "This page shows introduction to the project ",
                    "CoreMet.",
                    html.Br(),"For detailed information of the procedure, see ",
                    html.A("here", href="https://academic.oup.com/bib/article/24/4/bbad189/7176311", target="_blank")
                ]),
                id=idFunc("moreInfoCollapse"),
                is_open=False,
            ),
            dbc.Tooltip("More info.", target=idFunc("moreInfoIcon"), className='ms-1')
            ],
            fluid=True,
            className="py-1 bg-light rounded-3",
        ),
        className="p-0 my-1",
    )
    return header

def make_MPIDBHeader(idFunc):
    """
    Makes the header for the PV page
    """
    # Main Title
    header = html.Div(
        dbc.Container([
            html.Div([
                html.H2("Enzymatic reaction link database", className="display-4"),
                dbc.Button("Cite", id=idFunc("btn_citeHeader"),color="success", outline=True)
            ], className='d-flex justify-content-between align-items-center mb-0'),
            
            html.Hr(className="mt-0 mb-1"),
            html.Div([
                html.P("An atlas for Metabolite-Protein Interaction across ten organisms"),
                html.H4(id=idFunc('moreInfoIcon'), className="fa-solid fa-book-open ms-3 mt-1 primary")
            ], className='d-flex mb-0'),
            dbc.Collapse(
                html.Div([
                    "This page shows comprehensive information of metabolite-protein interaction information in ten different ",
                    "organisms", html.Br(), "Data was collected come from HMDB, SMPDB, ChEMBL, KEGG databases.", 
                    html.Br(),"For detailed information of the procedure, see ",
                    html.A("here", href="https://academic.oup.com/bib/article/24/4/bbad189/7176311", target="_blank")
                ]),
                id=idFunc("moreInfoCollapse"),
                is_open=False,
            ),
            dbc.Tooltip("More info.", target=idFunc("moreInfoIcon"), className='ms-1')
            ],
            fluid=True,
            className="py-1 bg-light rounded-3",
        ),
        className="p-0 my-1",
    )
    return header

def make_DockingHeader(idFunc):
    """
    Makes the header for the Genes page
    """
    # Main Title
    header = html.Div(
        dbc.Container([
            html.Div([
                html.H2("Visualization of Interaction Prediction by Molecular Docking", className="display-4"),
                dbc.Button("Cite", id=idFunc("btn_citeHeader"),color="success", outline=True)
            ], className='d-flex justify-content-between align-items-center mb-0'),
            
            html.Hr(className="mt-0 mb-1"),
            html.Div([
                html.P("Explore interaction results by molecular docking"),
                html.H4(id=idFunc('moreInfoIcon'), className="fa-solid fa-book-open ms-3 mt-1 primary")
            ], className='d-flex mb-0'),
            dbc.Collapse(
                html.Div([
                    "This page shows interactive visualizations of interaction results.",
                    html.Br(), 
                    html.Br(),"For detailed information of the procedure, see ",
                    html.A("here", href="https://academic.oup.com/bib/article/24/4/bbad189/7176311", target="_blank")
                ]),
                id=idFunc("moreInfoCollapse"),
                is_open=False,
            ),
            dbc.Tooltip("More info.", target=idFunc("moreInfoIcon"), className='ms-1')
            ],
            fluid=True,
            className="py-1 bg-light rounded-3",
        ),
        className="p-0 my-1",
    )
    return header

def make_introductionText():
    """
    Makes the introduction text for the introduction page
    """
    text = html.Div([
        html.P("Identification of metabolite-protein interactions (MPIs) is essential for understanding enzymatic reactions and metabolic networks but remains a significant challenge due to the complexity of these interactions. In this study we have developed CoreMet, a web server that utilizes GraphSAGE-based graph neural networks to predict MPIs and additional interaction types including metabolite–disease, metabolite–microbe, and metabolite–drug interactions. The CoreMet server provides a user-friendly interface that enables rapid prediction and visualization of MPIs and molecular docking results across a wide range of species and metabolites, including nucleic acids, lipids, organic acids, carbohydrates, and drugs. This code-free web server allows researchers, even those without expertise in computer programming or bioinformatics, to predict enzymatic reaction links efficiently and explore novel metabolic pathways."
               ),
        dbc.Card(
        [
            dbc.CardImg(src="/assets/workflow.png", top=True),
            dbc.CardBody(
                html.P("Workflow of CoreMet Framework", className="card-text")
            ),
        ],
        style={"txt_align": "justify"},
        ),
    ],        
    style={"txt_align": "justify",
         },)
    return text

import pandas as pd

def make_DatabaseInfo(idFunc,df):
    database = html.Div([
        html.H6(" ", className='my-1'),
        dash_table.DataTable(
        df.to_dict('records'),
        [{"name": i, "id": i} for i in df.columns],
        sort_action="native",
        sort_mode="multi",
        filter_action="native",
        filter_options={"placeholder_text": "Filter column..."},
        page_size=20,
        style_cell={
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
        'maxWidth': 0,
        'font-family': "Arial, Helvetica, sans-serif",
        'textAlign': 'left',
        },
        tooltip_data=[
        {
            column: {'value': str(value), 'type': 'markdown'}
            for column, value in row.items()
        } for row in df.to_dict('records')
        ],
        tooltip_duration=2000,
        style_table={'overflowX': 'auto',
                     'font-family': "Arial, Helvetica, sans-serif"},)  
    ])
    
    return database


def make_InteractionHeader(idFunc):
    """
    Makes the header for the WFA page
    """
    # Main Title
    header = html.Div(
        dbc.Container([
            html.Div([
                html.H2("Metabolite-Protein Interaction Prediction", className="display-4"),
                dbc.Button("Cite", id=idFunc("btn_citeHeader"),color="success", outline=True)
                ], className='d-flex justify-content-between align-items-center mb-0'
            ),
            html.Hr(className="mt-0 mb-1"),
            html.Div([
                html.P("Explore the relationship between metabolites and proteins"),
                html.H4(id=idFunc('moreInfoIcon'), className="fa-solid fa-book-open ms-3 mt-1 primary")
            ], className='d-flex mb-0'),
            dbc.Collapse(
                html.Div([
                    "This page shows the key process for predicting metabolite-protein interaction, ",
                    "including the upload of metabolite and protein data, the prediction of interaction, ",
                    "and the visualization of the result.",
                    html.Br(), 
                    html.Br(),"For detailed information of the procedure, see ",
                    html.A("here", href="https://academic.oup.com/bib/article/24/4/bbad189/7176311", target="_blank")
                ]),
                id=idFunc("moreInfoCollapse"),
                is_open=False,
            ),
            dbc.Tooltip("More info.", target=idFunc("moreInfoIcon"), className='ms-1')
            ],
            fluid=True,
            className="py-1 bg-light rounded-3",
        ),
        className="p-0 my-1",
    )
    return header

def make_AreasChecklist(idFunc,coarseDict):
    checklist = dbc.Checklist(
        id=idFunc('chklst_areasCorr'),
        options=coarseDict,
        value=[x['value'] for x in coarseDict],
        inline=True,
        className='ms-5'
    )
    return checklist

def make_InteractionUploadMenu(idFunc):
    """
    Makes the layout for the upload menu in the interaction page
    """
    uploadMenu = html.Div([
        html.H6("Upload the metabolite data"),
        dcc.Upload(
            id=idFunc('upload-metabolite-data'),
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Files')
            ]),
            style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            # Allow multiple files to be uploaded
            multiple=True
        ),
        html.H6("Or enter metabolite data manually:"),
        dcc.Textarea(
            id=idFunc('metabolite-prefill-textarea'),
            value='1-methyl-L-histidine,'
            'HMDB0000001,'
            'CN1C=NC(C[C@H](N)C(O)=O)=C1',
            style={
                'width': '100%',
                'height': '100px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'left',
                'margin': '10px'
            }
        ),
        html.Div(id=idFunc('metabolite-file')),
        html.H6("   "), 
        html.H6("Upload the protein data"),
        dcc.Upload(
            id=idFunc('upload-protein-data'),
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Files')
            ]),
            style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            # Allow multiple files to be uploaded
            multiple=True
        ),
        html.Div(id=idFunc('output-protein-data-upload')),
        html.H6("Or enter protein data manually:"),
        dcc.Textarea(
            id=idFunc('protein-prefill-textarea'),
            value='A0A075B6H7,'
            'Probable non-functional immunoglobulin kappa variable 3-7,'
            'IGKV3-7,'
            'Homo sapiens (Human),'
            'MEAPAQLLFLLLLWLPDTTREIVMTQSPPTLSLSPGERVTLSCRASQSVSSSY'
            'LTWYQQKPGQAPRLLIYGASTRATSIPARFSGSGSGTDFTLTISSLQPEDFAVYYCQQDYNLP',
            style={
                'width': '100%',
                'height': '100px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'left',
                'margin': '10px'
            }
        ),
        html.Div(id=idFunc('protein-file')),
        dcc.Store(id=idFunc('store-metabolite-data-upload'),storage_type='local',data={}),
        dcc.Store(id=idFunc('store-protein-data-upload'),storage_type='local',data={}),
        dbc.Button('Submit', id=idFunc('upload-button'), n_clicks=0),
    ])
    return uploadMenu

def make_InteractionTable(idFunc):
    uploadTable = html.Div([
        html.Div(id=idFunc('table-metabolite-data-upload')),
        html.Div(id=idFunc('table-protein-data-upload')),
    ])
    return uploadTable

def make_MPISelectionMenu(idFunc, genome_dict):
    """
    Makes the left-side menu with dropdowns for the histogram of the multiple staining metrics
    """
    menu = html.Div([
        html.H6(["Select a genome:"],className='my-1'),
        html.Div([
            html.Div([
                dcc.Dropdown(
                    id=idFunc('drpD_genome_Select'),
                    options=list(genome_dict.keys()),
                    value='Homo sapiens', # Default ID for Aggrecan
                    multi=False,
                    clearable=False,
                )],
                style={'flex-grow':'1'},
            ),

            html.Div([
                dbc.Button([html.I(className="fa-solid fa-info")],id=idFunc('btn_info'))], 
            )],
            className='mt-1 mb-1', 
            style={'display':'flex','flex-direction':'row', 'column-gap':'5px'}
        ),

        html.Div([
        dbc.Button("Open Organism Genome Info",
            id=idFunc("btn_openTabDiffuse"),
            className="mb-1",
            size='sm',
            outline=True,
            color="primary",
        ),
        dbc.Collapse(
            id=idFunc("collps_Tab"),
            is_open=False,
        )],className='mt-3 mb-5'),

        html.H6(["Enable protein similarity match?"],className='my-1'),
        
        dcc.Dropdown(
            id=idFunc('drpD_matchSelector'),
            options = ['Yes','No'],
            value='Yes', # Defaults to Isocortex
            multi = False,
            clearable=False,
            className='my-1 mb-3'
        ),

        dbc.Button('Submit', id=idFunc('feature-button'), n_clicks=0),
        dcc.Store(id=idFunc('store-molecular-feature'),storage_type='local',data={}),
        dcc.Store(id=idFunc('store-mpi-result-table'),storage_type='local',data={}),
        # TOOLTIPS
        dbc.Tooltip("Select the genome that you want to perform metabolite-protein interaction prediction.",
            target=idFunc("drpD_genome_Select"),placement="right"
        ),
        dbc.Tooltip("Select whether you want perform protein similarity search if target protein was not available.",target=idFunc("drpD_matchSelector"),placement="right"),
    ])
    return menu

def MPI_CollapsableTable(idFunc):
    """
    Makes the collapsable tabular data section
    """
    collapsTable = html.Div([
        dbc.Button("Show Interaction Prediction Result",
            id=idFunc("btn_openTabMPI"),
            className="mb-1",
            color="primary",
        ),
        dbc.Collapse(
            id=idFunc("collps_MPI_Tab"),
            is_open=False,
        )],
        className='mt-3'
    )
    return collapsTable

def make_MPIVisualizationMenu(idFunc, mets_dict, protein_dict):
    """
    Makes the left-side menu with dropdowns for the histogram of the multiple staining metrics
    """
    menu = html.Div([
        html.H6(["Select a Metabolite ID:"],className='my-1'),
        dcc.Dropdown(
            id=idFunc('drpD_metSelect'),
            options=list(mets_dict.keys()),
            value='HMDB0001488', # Default ID
            multi=False,
            clearable=False,
        ),
        html.H6(["\n\n\n"],className='my-1'),
        html.H6(["Select a Protein:"],className='my-1'),
        dcc.Dropdown(
            id=idFunc('drpD_ProteinSelect'),
            options = list(protein_dict.keys()),
            value='P53184', # Defaults
            multi = False,
            clearable=False,
            className='my-1 mb-3'
        ),
        dbc.Button('Submit', id=idFunc('dock-button'), n_clicks=0),
        dcc.Store(id=idFunc('pdb'),storage_type='local',data={}),
        dcc.Store(id=idFunc('poses'),storage_type='local',data={}),
        # TOOLTIPS
        dbc.Tooltip("Select the metabolite that interact protein with.",
            target=idFunc("drpD_metSelect"),placement="right"
        ),
        dbc.Tooltip("Select the protein to use for interaction visualization by molecular docking.",target=idFunc("drpD_ProteinSelect"),placement="right"),
    ])
    return menu


def make_InteractionSelectionMenu(idFunc):
    """
    Makes the left-side menu of the correlation plot in the interaction page
    """

    # Create a list of dicts for the dropdown of WFA vs PV
    wfaPvDictList = [{'label':'WFA', 'value':'wfa'},{'label':'Parvalbumin', 'value':'pv'}]

    # Create a list of dicts for the dropdown of the different metrics
    metricsDictList = [
        {'label':'Cell Energy','value':'energy'},
        {'label':'Diffuse Fluorescence','value':'diffuseFluo'},
        {'label':'Cell Density','value':'density'},
        {'label':'Cell Intensity','value':'intensity'},
    ]

    menu = html.Div([

        # X AXIS
        dbc.Row([
            html.Div([
                html.H4("X axis", className='my-1 mt-2'),
                html.H6("Select a staining", className='mt-3 mb-1'),
                dcc.Dropdown(
                    id=idFunc('drpD_xStaining'),
                    options = wfaPvDictList,
                    value='wfa',
                    multi = False,
                    clearable=False,
                    className=''
                ),
                html.H6("Select a metric", className='mt-3 mb-1'),
                dcc.Dropdown(
                    id=idFunc('drpD_xMetric'),
                    options = metricsDictList,
                    value='energy',
                    multi = False,
                    clearable=False,
                    className='mb-2'
                ),

            ],className='my-3 border border-light rounded-3')
        ]),

        # Y AXIS
        dbc.Row([
            html.Div([
                html.H4("Y axis", className='my-1 mt-2'),
                html.H6("Select a staining", className='mt-3 mb-1'),
                dcc.Dropdown(
                    id=idFunc('drpD_yStaining'),
                    options = wfaPvDictList,
                    value='pv',
                    multi = False,
                    clearable=False,
                    className=''
                ),
                html.H6("Select a metric", className='mt-3 mb-1'),
                dcc.Dropdown(
                    id=idFunc('drpD_yMetric'),
                    options = metricsDictList,
                    value='energy',
                    multi = False,
                    clearable=False,
                    className='mb-2'
                ),

            ],className='my-3 border border-light rounded-3')
        ]),
        dbc.Row([
            html.Div([
            # dbc.Button("Reset", id="btn-reset", color="primary", className="my-1"),
            dbc.Switch(
                label="Z-Score",
                value=False,
                id=idFunc("switch_zScore"),
                # class_name="mt-1"
                )],
            className="d-flex"
            )
        ])
    ])
    return menu

def make_MetricsHistogramSelectionMenu(idFunc, coarseDict, midDict, fineDict, staining='wfa'):
    """
    Makes the left-side menu with dropdowns for the histogram of the multiple staining metrics
    """
    menu = html.Div([
        html.H6(["Select a metric to show:"],className='my-1'),

        html.Div([
            html.Div([
                dcc.Dropdown(
                    id=idFunc('drpD_histogMetric'),
                    options=getMetricsLabels(staining=staining),
                    value='energy',
                    multi = False,
                    clearable=False,
                )],
                style={'flex-grow':'1'},
            ),
            html.Div([
                dbc.Button([html.I(className="fa-solid fa-info")],id=idFunc('btn_info'))], 
            )],
            className='mt-1 mb-5', 
            style={'display':'flex','flex-direction':'row', 'column-gap':'5px'}
        ),
        html.H6(["Select a major brain subdivision:"],className='my-1'),
        dcc.Dropdown(
            id=idFunc('drpD_majorSubd'),
            options = coarseDict,
            value=315, # Defaults to Isocortex
            multi = False,
            className='my-1 mb-3'
        ),
        html.H6(["Add ",html.B("coarse-ontology")," region(s):"], className='my-1'),
        html.Div([
            html.Div([
                dcc.Dropdown(
                    id=idFunc('drpD_addCoarse'),
                    options = coarseDict,
                    multi = True,
                )],
                style={'flex-grow':'1'},
            ),
            html.Div([
                dbc.Button("All",id=idFunc('btn_allMajorDiff'))], 
            )],
            className='mt-1 mb-3', 
            style={'display':'flex','flex-direction':'row', 'column-gap':'5px'}
        ),
        html.H6(["Add ", html.B("mid-ontology"), " region(s):"],className='my-1'),
        dcc.Dropdown(
            id=idFunc('drpD_addMid'),
            options = midDict,
            multi = True,
            className="mt-1 mb-3"
        ),
        html.H6(["Add ", html.B("fine-ontology"), " region(s):"],className='my-1'),
        dcc.Dropdown(
            id=idFunc('drpD_addFine'),
            options = fineDict,
            multi = True,
            className="mt-1 mb-3"
        ),
        html.Div([
            # dbc.Button("Reset", id="btn-reset", color="primary", className="my-1"),
            dbc.Switch(
                label="Sort regions",
                value=False,
                id=idFunc("switch_sortDiff"),
                # class_name="mt-1"
                )],
            className="d-flex"
        ),

        # TOOLTIPS
        dbc.Tooltip("Add all the 12 major brain subdivisions.", target=idFunc("btn_allMajorDiff")),
        dbc.Tooltip(["Sort on intensity.", html.Br(), "Otherwise, sorted on region ontology."],
            target=idFunc("switch_sortDiff"),placement='bottom'
        ),
        dbc.Tooltip("Add all the regions of a selected major brain subdivision",
            target=idFunc("drpD_majorSubd"),placement="right"
        ),
        dbc.Tooltip("Major brain subdivisions",target=idFunc("drpD_addCoarse"),placement="right"),
        dbc.Tooltip("Regions at single-area level",target=idFunc("drpD_addMid"),placement="right"),
        dbc.Tooltip("Regions at single-layer level",target=idFunc("drpD_addFine"),placement="right"),
    ])
    return menu



def make_ColocalizationHistogramSelectionMenu(idFunc, coarseDict, midDict, fineDict):
    """
    Makes the left-side menu with dropdowns for the histogram of the multiple staining metrics
    """
    menu = html.Div([
        html.H6(["Select a colocalization metric:"],className='my-1'),

        html.Div([
            html.Div([
                dcc.Dropdown(
                    id=idFunc('drpD_Metric'),
                    options=[
                        {'label':'Percentage of PV-positive PNNs','value':'pvPositive_pnn'},
                        {'label':'Percentage of WFA-positive PV cells','value':'wfaPositive_pv'},
                    ],
                    value='pvPositive_pnn',
                    multi = False,
                    clearable=False,
                )],
                style={'flex-grow':'1'},
            ),
            html.Div([
                dbc.Button([html.I(className="fa-solid fa-info")],id=idFunc('btn_info'))], 
            )],
            className='mt-1 mb-5', 
            style={'display':'flex','flex-direction':'row', 'column-gap':'5px'}
        ),
        html.H6(["Select a major brain subdivision:"],className='my-1'),
        dcc.Dropdown(
            id=idFunc('drpD_majorSubd'),
            options = coarseDict,
            value=315, # Defaults to Isocortex
            multi = False,
            className='my-1 mb-3'
        ),
        html.H6(["Add ",html.B("coarse-ontology")," region(s):"], className='my-1'),
        html.Div([
            html.Div([
                dcc.Dropdown(
                    id=idFunc('drpD_addCoarse'),
                    options = coarseDict,
                    multi = True,
                )],
                style={'flex-grow':'1'},
            ),
            html.Div([
                dbc.Button("All",id=idFunc('btn_allMajorDiff'))], 
            )],
            className='mt-1 mb-3', 
            style={'display':'flex','flex-direction':'row', 'column-gap':'5px'}
        ),
        html.H6(["Add ", html.B("mid-ontology"), " region(s):"],className='my-1'),
        dcc.Dropdown(
            id=idFunc('drpD_addMid'),
            options = midDict,
            multi = True,
            className="mt-1 mb-3"
        ),
        html.H6(["Add ", html.B("fine-ontology"), " region(s):"],className='my-1'),
        dcc.Dropdown(
            id=idFunc('drpD_addFine'),
            options = fineDict,
            multi = True,
            className="mt-1 mb-3"
        ),
        html.Div([
            # dbc.Button("Reset", id="btn-reset", color="primary", className="my-1"),
            dbc.Switch(
                label="Sort regions",
                value=False,
                id=idFunc("switch_sortDiff"),
                # class_name="mt-1"
                )],
            className="d-flex"
        ),

        # TOOLTIPS
        dbc.Tooltip("Add all the 12 major brain subdivisions.", target=idFunc("btn_allMajorDiff")),
        dbc.Tooltip(["Sort in descending order.", html.Br(), "Otherwise, sorted on region ontology."],
            target=idFunc("switch_sortDiff"),placement='bottom'
        ),
        dbc.Tooltip("Add all the regions of a selected major brain subdivision",
            target=idFunc("drpD_majorSubd"),placement="right"
        ),
        dbc.Tooltip("Major brain subdivisions",target=idFunc("drpD_addCoarse"),placement="right"),
        dbc.Tooltip("Regions at single-area level",target=idFunc("drpD_addMid"),placement="right"),
        dbc.Tooltip("Regions at single-layer level",target=idFunc("drpD_addFine"),placement="right"),
    ])
    return menu

def make_GeneCorrSelectionMenu(idFunc, genesDf):
    """
    Makes the left-side menu with dropdowns for the histogram of the multiple staining metrics
    """
    menu = html.Div([
        html.H6(["Select a gene ID:"],className='my-1'),
        html.Div([
            html.Div([
                dcc.Dropdown(
                    id=idFunc('drpD_geneSelect'),
                    options=getGenesLabels(genesDf),
                    value=11382, # Default ID for Aggrecan
                    multi=False,
                    clearable=False,
                )],
                style={'flex-grow':'1'},
            ),

            html.Div([
                dbc.Button([html.I(className="fa-solid fa-info")],id=idFunc('btn_info'))], 
            )],
            className='mt-1 mb-1', 
            style={'display':'flex','flex-direction':'row', 'column-gap':'5px'}
        ),

        html.Div([
        dbc.Button("Open Gene Info",
            id=idFunc("btn_openTabDiffuse"),
            className="mb-1",
            size='sm',
            outline=True,
            color="primary",
        ),
        dbc.Collapse(
            id=idFunc("collps_Tab"),
            is_open=False,
        )],className='mt-3 mb-5'),

        html.H6(["Correlation with:"],className='my-1'),
        dcc.Dropdown(
            id=idFunc('drpD_metricSelector'),
            options = getMetricsForGenesLabels(),
            value='wfa_energy', # Defaults to Isocortex
            multi = False,
            clearable=False,
            className='my-1 mb-3'
        ),

        # TOOLTIPS
        dbc.Tooltip("Select the gene that you want to correlate staining metrics with.",
            target=idFunc("drpD_geneSelect"),placement="right"
        ),
        dbc.Tooltip("Select what staining metric to use for correlation.",target=idFunc("drpD_metricSelector"),placement="right"),
    ])
    return menu

def make_CollapsableTable(idFunc):
    """
    Makes the collapsable tabular data section
    """
    collapsTable = html.Div([
        dbc.Button("Open Tabular Data",
            id=idFunc("btn_openTabDiffuse"),
            className="mb-1",
            color="primary",
        ),
        dbc.Collapse(
            id=idFunc("collps_Tab"),
            is_open=False,
        )],
        className='mt-3'
    )
    return collapsTable

def make_AnatomicalExplorerSelectionMenu(idFunc, staining='wfa'):
    """
    Makes the layout for the left-side selection menu of the anatomical explorer
    """
    menu = html.Div([
        html.H6(["Select a metric to show:"],className='my-1'),
        html.Div([
            html.Div([
                dcc.Dropdown(
                    id=idFunc('drpD_anatomMetric'),
                    options= getMetricsLabels(staining=staining),
                    value='energy',
                    multi = False,
                    clearable=False,
                )],
                style={'flex-grow':'1'},
            ),
            html.Div([
                dbc.Button([html.I(className="fa-solid fa-info")],id=idFunc('btn_info_anat'))], 
            )],
            className='mt-1 mb-3', 
            style={'display':'flex','flex-direction':'row', 'column-gap':'5px'}
        ),

        html.H6(["Colormap:"],className='my-1'),
        dcc.Dropdown(
            id=idFunc('drpD_anatomCmap'),
            options = colormapDictListDropdown(),
            value='PuBu' if staining=='wfa' else 'Reds',
            multi = False,
            clearable=False,
            className="mt-1 mb-3"
        ),

        html.H6(["Antero-Posterior Axis:"],className='mt-5 mb-1'),
        dcc.Slider(0, 34, 1, value=10, id=idFunc('slider_ap'),
            marks={0:'Anterior',34:'Posterior'},
        ),

        # TOOLTIPS
        dbc.Tooltip("Visualization colormap.", target=idFunc("drpD_anatomCmap")),
        dbc.Tooltip("Set minimum and maximum values.", target=idFunc("slider_clims")),
    ])
    return menu

def make_CitationOffCanvas(idFunc):
    """
    Makes the layout for the offcanvas menu for citations
    """
    offcanvas = dbc.Offcanvas(
        html.P([
            "If you found this web server to be useful in your research, please consider ", 
            html.B("citing us."),
            html.Br(),
            dbc.Card(
            [
                dbc.CardImg(
                    src="/assets/network.jpg",
                    top=True,
                    style={"opacity": 0.4},
                ),
                dbc.CardImgOverlay(
                    dbc.CardBody(
                        [
                            html.H4("Briefings in bioinformatics", className="card-title"),
                            html.P(
                                "Check out our method paper here!",
                                className="card-text",
                            ),
                            dcc.Link([dbc.Button("Go to Paper", color="primary")],
                                href="https://academic.oup.com/bib/article/24/4/bbad189/7176311",
                                target="_blank"),
                            
                        ],
                    ),
                ),
            ],
            style={"width": "18rem"},className='mt-5'),
            # License banner removed to avoid duplicate IDs - it appears at the bottom of the page
            html.Div([
                html.Hr(className="mt-2 mb-2"),
                html.P(["Web-app developed by ",
                    dcc.Link("Cheng Wang",
                        href="https://scholar.google.com/citations?user=UAZhchQAAAAJ&hl=en",
                        target="_blank", className="me-3"),
                "See source code ",
                    dcc.Link("here",
                        href="https://github.com/mmetalab/mpi-vgae-web",
                        target="_blank"),
                ]),
                html.A([
                    html.Img([], alt="Creative Commons License",  
                        src="https://i.creativecommons.org/l/by/4.0/88x31.png")], 
                    rel="license", href="http://creativecommons.org/licenses/by/4.0/", className="border-width:0 me-2"),
                "This work is licensed under a ",
                html.A(["Creative Commons Attribution 4.0 International License"], rel='license', href="http://creativecommons.org/licenses/by/4.0/")
            ], className='pt-5'),  # No ID to avoid duplication
        ]),
        
        id=idFunc("offCanv_cite"),
        title="Cite us",
        is_open=False,
    )
    return offcanvas

def make_AboutUsOffCanvas(idFunc):
    offcanvas = dbc.Offcanvas(
        html.P([
            "This work was produced in Cheng Wang's Lab thanks to the support of the following fundings: ",
            html.Br(),
            dbc.ListGroup([
                dbc.ListGroupItem(
                    html.Div([
                        "National Natural Science Foundation of China",
                        html.A(className="fa-solid fa-arrow-up-right-from-square px-3", href='https://www.nsfc.gov.cn/english/site_1/index.html', target='_blank')
                    ], className='d-flex w-100 justify-content-between align-items-center'),
                    className='rounded-3 py2'),
                dbc.ListGroupItem(
                    html.Div([
                        "Shandong Natural Science Foundation",
                        html.A(className="fa-solid fa-arrow-up-right-from-square px-3", href='http://kjt.shandong.gov.cn/', target='_blank')
                    ], className='d-flex w-100 justify-content-between align-items-center'),
                    className='rounded-3 py2'),
            ], flush=True, className='mt-5'),
        ]),
        
        id=idFunc("offCanv_abtUs"),
        title="More about us",
        is_open=False,
        placement='end',
    )
    return offcanvas

def make_MetricInfoModal(idFunc):

    diffuseCard = dbc.Card([
        dbc.CardBody([
            html.H5("Diffuse Fluorescence (A.U.)", className="card-title primary"),
            html.P(html.I("'How much staining is there?'")),
            html.Hr(className="mb-1"),
            html.P(["For each mouse, we calculate the ", html.B("average fluorescence"), 
                " intensity across all the pixels belonging to a brain region. " + 
                "Then, we normalize this value dividing it by the average fluorescence intensity " + 
                "of the entire brain."],
                className="card-text")]
        )],
        style={'height':'100%'},
        outline=True,
        color='info',
    )

    densityCard = dbc.Card([
        dbc.CardBody([
            html.H5(["Cell Density (cells/mm", html.Sup(2) ,")"], className="card-title"),
            html.P(html.I("'How many cells/PNNs are there?'")),
            html.Hr(className="mb-1"),
            html.P(["For each mouse, we calculate the ", html.B(["number of cells per mm",html.Sup(2)]), 
                " In a brain region. " + 
                "This value is the cell density."],
                className="card-text")]
        )],
        style={'height':'100%'},
        outline=True,
        color='info',
    )

    intensityCard = dbc.Card([
        dbc.CardBody([
            html.H5("Cell Intensity (A.U.)", className="card-title"),
            html.P(html.I("'How intense are the cells/PNNs here?'")),
            html.Hr(className="mb-1"),
            html.P(["For each cell, we measure the average pixel intesity value (", html.I(["range: 0-1"]), ") of ",
                "all pixels belonging to that cell.", html.Br(), "For each brain region, we compute the ",
                "average of the intensity of all the cells in that region."],
                className="card-text")],        
        )],
        style={'height':'100%'},
        outline=True,
        color='info'
    )

    energyCard = dbc.Card([
        dbc.CardBody([
            html.H5("Cell Energy (A.U.)", className="card-title"),
            html.P(html.I("'What's the overall strenght of cells/PNNs here?'")),
            html.Hr(className="mb-1"),
            html.P(["Can be thought of as a measure of cell density, weighted by intensity. ",html.Br(),
                "For each region, energy is defined as the sum of the cell intensity ",
                "of all the cells in that region, divided by the total surface area ",
                "Then, we normalize this value dividing it by the energy of the entire brain."],
                className="card-text")]
        )],
        style={'height':'100%'},
        outline=True,
        color='info',
    )

    modal = dbc.Modal([
            dbc.ModalHeader(dbc.ModalTitle("Metrics analyzed")),
            dbc.ModalBody([
                dbc.Row([
                    dbc.Col([diffuseCard], width=6),
                    dbc.Col([densityCard]),
                ], className=''),
                dbc.Row([
                    dbc.Col([intensityCard], width=6),
                    dbc.Col([energyCard]),
                ], className='my-3')]
        )],
        size='lg',
        id=idFunc('modal_info')
    )

    return modal

def make_ColocInfoModal(idFunc):

    pvPositivePnnCard = dbc.Card([
        dbc.CardBody([
            html.H5("PV-positive PNNs", className="card-title primary"),
            html.P(html.I("'What's the percentage of PNNs that surround a PV cell?'")),
            html.Hr(className="mb-1"),
            html.P(["For each mouse, we calculate the ", html.B("average fluorescence"), 
                " intensity across all the pixels belonging to a brain region. " + 
                "Then, we normalize this value dividing it by the average fluorescence intensity " + 
                "of the entire brain."],
                className="card-text")]
        )],
        style={'height':'100%'},
        outline=True,
        color='info',
    )

    wfaPositivePvCard = dbc.Card([
        dbc.CardBody([
            html.H5(["Cell Density (cells/mm", html.Sup(2) ,")"], className="card-title"),
            html.P(html.I("'What's the percentage of PV cells that are surrounded by a PNN?''")),
            html.Hr(className="mb-1"),
            html.P(["For each mouse, we calculate the ", html.B(["number of cells per mm",html.Sup(2)]), 
                " In a brain region. " + 
                "This value is the cell density."],
                className="card-text")]
        )],
        style={'height':'100%'},
        outline=True,
        color='info',
    )

    modal = dbc.Modal([
            dbc.ModalHeader(dbc.ModalTitle("Colocalization metrics analyzed")),
            dbc.ModalBody([
                dbc.Row([
                    dbc.Col([pvPositivePnnCard], width=6),
                    dbc.Col([wfaPositivePvCard]),
                ], className='')]
        )],
        size='lg',
        id=idFunc('modal_info')
    )

    return modal

def make_GeneInfoModal(idFunc):

    genesCard = dbc.Card([
        dbc.CardBody([
            # html.H5("Gene ID", className="card-title primary"),
            html.P(html.I("Genome used for the MPI prediction")),
            html.Hr(className="mb-1"),
            html.P(["As a prediction of metabolite-protein interaction for ", html.B("organism specific"), 
                " as defined in Wang et al. 2023. "],
                className="card-text")]
        )],
        style={'height':'100%'},
        outline=True,
        color='info',
    )

    modal = dbc.Modal([
            dbc.ModalHeader(dbc.ModalTitle("Select organism genome")),
            dbc.ModalBody([
                dbc.Row([
                    dbc.Col([genesCard]),
                ], className='')]
        )],
        size='lg',
        id=idFunc('modal_info')
    )

    return modal

def make_CC_licenseBanner(idFunc):
    banner = []

    banner = html.Div([
        html.Hr(className="mt-2 mb-2"),
        html.P(["Web-app developed by ",
            dcc.Link("Cheng Wang",
                href="https://scholar.google.com/citations?user=UAZhchQAAAAJ&hl=en",
                target="_blank", className="me-3"),
        "See source code ",
            dcc.Link("here",
                href="https://github.com/mmetalab/mpi-vgae-web",
                target="_blank"),
        ]),
        html.A([
            html.Img([], alt="Creative Commons License",  
                src="https://i.creativecommons.org/l/by/4.0/88x31.png")], 
            rel="license", href="http://creativecommons.org/licenses/by/4.0/", className="border-width:0 me-2"),
        "This work is licensed under a ",
        html.A(["Creative Commons Attribution 4.0 International License"], rel='license', href="http://creativecommons.org/licenses/by/4.0/")
    ], id=idFunc("licenseBanner"), className='pt-5')
    return banner

# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------

def colormapDictListDropdown():
    """
    Creates a list of dicts to fill a dropdown to select different colormaps
    """
    cmapDictList = [
        dict(label='Blue',value='PuBu'),
        dict(label='Red',value='Reds'),
        dict(label='Green',value='Greens'),
        dict(label='Gray',value='Greys'),
        dict(label='Viridis',value='viridis'),
        dict(label='Magma',value='magma'),
    ]
    return cmapDictList

def getMetricsLabels(staining='wfa'):
    if staining == 'wfa':
        labels = [
            {'label':'PNN Energy','value':'energy'},
            {'label':'Diffuse Fluorescence','value':'diffuseFluo'},
            {'label':'PNN Density','value':'density'},
            {'label':'PNN Intensity','value':'intensity'},
        ]
    elif staining == 'pv':
        labels = [
            {'label':'PV Energy','value':'energy'},
            {'label':'Diffuse Fluorescence','value':'diffuseFluo'},
            {'label':'PV Density','value':'density'},
            {'label':'PV Intensity','value':'intensity'},
        ]

    return labels

def getGenesLabels(genesDf):
    genesDf = genesDf.sort_values(by='gene_acronym')
    labels = []
    for idx, row in genesDf.iterrows():
        temp = {'label':idx,'value':row['gene_AGEA_id']}
        labels.append(temp)
    return labels

def getMetricsForGenesLabels():
    labels = [
        {'label':'PNN Energy','value':'wfa_energy'},
        {'label':'WFA Diffuse Fluorescnece','value':'wfa_diffuseFluo'},
        {'label':'PV Energy','value':'pv_energy'},
    ]
    return labels