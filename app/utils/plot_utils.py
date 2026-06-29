"""
Plotting utility functions for CoreMet Application
"""

import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from typing import Dict, Any, Optional


def create_scatter_plot(data: pd.DataFrame, 
                       x_col: str = 'UMAP-1', 
                       y_col: str = 'UMAP-2',
                       color_col: str = 'Molecular type',
                       title: str = "Molecular Feature Visualization") -> go.Figure:
    """
    Create a scatter plot for molecular features
    
    Args:
        data: DataFrame with plot data
        x_col: Column name for x-axis
        y_col: Column name for y-axis
        color_col: Column name for color mapping
        title: Plot title
        
    Returns:
        Plotly figure object
    """
    if data.empty:
        return create_empty_plot(title)
    
    fig = px.scatter(
        data, 
        x=x_col, 
        y=y_col,
        color=color_col,
        hover_data=['Name'] if 'Name' in data.columns else None,
        title=title
    )
    
    # Customize layout
    fig.update_layout(
        title={
            'text': title,
            'y': 0.9,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        font=dict(family="Arial, Helvetica, sans-serif"),
        height=600,
        showlegend=True
    )
    
    # Customize axes
    fig.update_xaxes(title=x_col)
    fig.update_yaxes(title=y_col)
    
    return fig


def update_scatter_plot(fig: go.Figure, 
                       data: pd.DataFrame,
                       x_col: str = 'UMAP-1',
                       y_col: str = 'UMAP-2',
                       color_col: str = 'Molecular type') -> go.Figure:
    """
    Update existing scatter plot with new data
    
    Args:
        fig: Existing plotly figure
        data: New data to plot
        x_col: Column name for x-axis
        y_col: Column name for y-axis
        color_col: Column name for color mapping
        
    Returns:
        Updated plotly figure
    """
    if data.empty:
        return create_empty_plot("No data available")
    
    # Clear existing data
    fig.data = []
    
    # Add new data
    fig = px.scatter(
        data,
        x=x_col,
        y=y_col,
        color=color_col,
        hover_data=['Name'] if 'Name' in data.columns else None
    )
    
    return fig


def create_empty_plot(title: str = "No Data Available") -> go.Figure:
    """
    Create an empty plot with a message
    
    Args:
        title: Title for the empty plot
        
    Returns:
        Empty plotly figure
    """
    fig = go.Figure()
    
    fig.add_annotation(
        text="No data available for visualization",
        xref="paper",
        yref="paper",
        x=0.5,
        y=0.5,
        showarrow=False,
        font=dict(size=16, color="gray")
    )
    
    fig.update_layout(
        title=title,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        height=400
    )
    
    return fig


def create_prediction_results_table(data: pd.DataFrame) -> Dict[str, Any]:
    """
    Create table configuration for prediction results
    
    Args:
        data: DataFrame with prediction results
        
    Returns:
        Dictionary with table configuration
    """
    if data.empty:
        return {
            'data': [],
            'columns': [],
            'style': {}
        }
    
    # Define column configurations
    columns = []
    for col in data.columns:
        column_config = {
            'id': col,
            'name': col,
            'type': 'text'
        }
        
        # Special formatting for specific columns
        if 'Prediction Score' in col:
            column_config['type'] = 'numeric'
            column_config['format'] = {'specifier': '.5f'}
        elif 'Existing' in col:
            column_config['type'] = 'text'
        
        columns.append(column_config)
    
    # Style configuration
    style = {
        'textAlign': 'left',
        'minWidth': '150px',
        'width': '150px',
        'maxWidth': '200px',
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
    }
    
    return {
        'data': data.to_dict('records'),
        'columns': columns,
        'style': style,
        'style_table': {
            'minWidth': '100%',
            'height': 400,
            'overflowX': 'auto'
        }
    }


def create_database_info_table(data: pd.DataFrame) -> Dict[str, Any]:
    """
    Create table configuration for database information
    
    Args:
        data: DataFrame with database information
        
    Returns:
        Dictionary with table configuration
    """
    if data.empty:
        return {
            'data': [],
            'columns': [],
            'style': {}
        }
    
    columns = [{"name": i, "id": i} for i in data.columns]
    
    style = {
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
        'maxWidth': 0,
        'font-family': "Arial, Helvetica, sans-serif",
        'textAlign': 'left',
    }
    
    return {
        'data': data.to_dict('records'),
        'columns': columns,
        'style': style,
        'style_table': {
            'overflowX': 'auto',
            'font-family': "Arial, Helvetica, sans-serif"
        },
        'tooltip_data': [
            {
                column: {'value': str(value), 'type': 'markdown'}
                for column, value in row.items()
            } for row in data.to_dict('records')
        ],
        'tooltip_duration': 2000
    }


def create_upload_preview_table(data: pd.DataFrame, data_type: str) -> Dict[str, Any]:
    """
    Create table configuration for upload preview
    
    Args:
        data: DataFrame with uploaded data
        data_type: Type of data ('metabolite' or 'protein')
        
    Returns:
        Dictionary with table configuration
    """
    if data.empty:
        return {
            'data': [],
            'columns': [],
            'style': {}
        }
    
    columns = [{"name": i, "id": i} for i in data.columns]
    
    # Adjust column widths based on data type
    if data_type == 'metabolite':
        column_widths = {
            'Metabolite Name': '200px',
            'HMDB ID': '120px',
            'SMILES': '300px'
        }
    else:  # protein
        column_widths = {
            'UniprotID': '120px',
            'Protein Name': '200px',
            'Gene Name': '120px',
            'Organism': '150px',
            'Sequence': '300px'
        }
    
    style = {
        'textAlign': 'left',
        'minWidth': '100px',
        'width': '100px',
        'maxWidth': '300px',
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
    }
    
    return {
        'data': data.to_dict('records'),
        'columns': columns,
        'style': style,
        'style_table': {
            'minWidth': '100%',
            'height': 270
        }
    }
