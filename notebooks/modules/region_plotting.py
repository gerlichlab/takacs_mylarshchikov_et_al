from dataclasses import dataclass
from typing import Any, Callable, Literal, Optional

import bbi
import bioframe as bf
import cooler
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colors as mcolors
from matplotlib import transforms as mtransforms
from matplotlib.patches import Patch, Rectangle
from numpy.typing import ArrayLike


@dataclass
class MtxTrack:
    data: ArrayLike
    vmin: float
    vmax: float
    sep: int
    norm: Any
    cmap: Any
    trim: bool = True


@dataclass
class LineTrack:
    data: ArrayLike
    vmin: float
    vmax: float
    color: str


@dataclass
class FillTrack:
    data: ArrayLike
    vmin: float
    vmax: float
    neg_color: str
    pos_color: str


@dataclass
class IntervalTrack:
    data: pd.DataFrame
    color: str
    edgecolor: str
    pad: int = 0


@dataclass
class MultiIntervalTrack:
    data: pd.DataFrame
    colors: dict[str, str]
    edgecolors: dict[str, str]
    hue: str
    stratify: bool = False
    pad: int = 0


@dataclass
class TADTrack:
    data: pd.DataFrame
    facecolor: Any
    edgecolor: Any
    facealpha: float
    edgealpha: float
    orient: Literal['up', 'down']

    def __post_init__(self):
        orient_options = {'up', 'down'}
        if self.orient not in orient_options:
            raise ValueError(f"Invalid orient value: {self.orient}. Must be one of {orient_options}")
        

@dataclass
class RulerTrack:
    loc: Literal['top', 'bottom']
    locator: Any
    formatter: Any

    def __post_init__(self):
        loc_options = {'top', 'bottom'}
        if self.loc not in loc_options:
            raise ValueError(f"Invalid loc value: {self.loc}. Must be one of {loc_options}")


class RegionPlotter:
    """
    RegionPlotter is a class for visualizing genomic regions with various types of tracks, 
    including Hi-C heatmaps, signal tracks, and interval tracks. It provides methods to 
    add and customize tracks, as well as generate plots for the specified genomic region.
    Attributes:
        chrom (str): The chromosome identifier for the genomic region.
        start (int): The start position of the genomic region.
        end (int): The end position of the genomic region.
        tracks (list): A list of tracks added to the plot.
        heights (list): A list of heights corresponding to the tracks.
    Methods:
        __init__(chrom, start, end):
            Initializes the RegionPlotter with the specified chromosome, start, and end positions.
        _read_mtx(clr, pad=0):
            Reads a Hi-C matrix from a Cooler object for the specified genomic region with optional padding.
        _add_mtx(mtx, vmin, vmax, sep, norm, cmap, height=5, trim=False):
            Adds a Hi-C matrix track to the plot with specified visualization parameters.
        add_hic(clr_path, resolution, vmin, vmax, sep, cmap, norm='log', height=5):
            Adds Hi-C data from a Cooler file to the plot.
        _read_array(signal_path, binsize, mean_subtract=False):
            Reads a signal array from a BigWig file for the specified genomic region.
        _add_array(array, kind, height=2, **kwargs):
            Adds a signal array to the plot as a line or filled track.
        add_line(signal_path, binsize, vmin=None, vmax=None, color='k', mean_subtract=False, height=2):
            Adds a signal as a line track to the plot.
        add_fill(signal_path, binsize, vmin=None, vmax=None, neg_color='grey', pos_color='k', mean_subtract=False, height=2):
            Adds a signal as a filled track to the plot.
        _read_bed(bed_path):
            Reads a BED-compatible file and returns its contents as a DataFrame.
        add_df(self, df, color='k', edgecolor=None, pad=0, height=1):
            Adds a DataFrame with genomic intervals to the plot.
        add_df_multi(self, df, hu, colors, edgecolors=None, pad=0, stratify=False, height=1):
            Adds a DataFrame with genomic intervals and multiple categories to the plot.
        add_bed(self, bed_path, color='k', edgecolor=None, pad=0, height=1):
            Adds a BED-compatible file with genomic intervals to the plot.
        _plot_hic(track, ax, cax):
            Plots a Hi-C heatmap on the given axes.
        _plot_line(track, ax):
            Plots a line track on the given axis.
        _plot_fill(track, ax):
            Plots a filled track on the given axis.
        _plot_df(track, ax):
            Plots interval data from an IntervalTrack on the given axis.
        _plot_df_multi(track, ax, lax):
            Plots interval data with multiple categories from a MultiIntervalTrack on the given axes.
        plot(figsize, dpi, width_ratios=(40, 1), tick_all=True, tick_step=1_000_000):
            Generates a plot for the tracks with specified layout and formatting.
    """

    def __init__(self, chrom: str, start: int, end: int):
        """
        Initializes a region with specified chromosome, start, and end positions.
        Args:
            chrom (str): The chromosome identifier.
            start (int): The start position of the region.
            end (int): The end position of the region.
        Attributes:
            chrom (str): The chromosome identifier.
            start (int): The start position of the region.
            end (int): The end position of the region.
            tracks (list): A list to store associated tracks for the region.
            heights (list): A list to store heights corresponding to tracks.
        """
        self.chrom: str = chrom
        self.start: int = start
        self.end: int = end
        self.tracks = list()
        self.heights = list()

    def _read_mtx(self, clr: cooler.Cooler, pad: int = 0) -> ArrayLike:
        """
        Reads a matrix from a Cooler object for a specified genomic region with optional padding.
        Args:
            clr (cooler.Cooler): A Cooler object representing the Hi-C contact matrix.
            pad (int, optional): The number of base pairs to pad around the specified region. 
                Defaults to 0.
        Returns:
            numpy.ndarray: A 2D array representing the extracted matrix for the specified region.
        """
        mtx = clr.matrix().fetch((self.chrom, self.start - pad, self.end + pad))
        return mtx

    def _add_mtx(self, mtx: ArrayLike, vmin: float, vmax: float, sep: int, cmap: Any, norm: Any, height: int = 5, trim: bool = False):
        """
        Adds a matrix track to the plotting object.
        Args:
            mtx (ArrayLike): The matrix data to be added as a track.
            vmin (float): The minimum value for the colorbar.
            vmax (float): The maximum value for the colorbar.
            sep (int): The maximum genomic separation (in bins) to display on the image.
            cmap (Any): The colormap to be used for the track.
            norm (Any): The normalization function or object to display the matrix.
            height (int, optional): The height of the track. Defaults to 5.
            trim (bool, optional): Whether to trim the matrix data for plotting. Defaults to False.
        """
        track = MtxTrack(mtx, vmin, vmax, sep=sep, norm=norm, cmap=cmap, trim=trim)
        self.tracks.append(track)
        self.heights.append(height)

    def add_hic(self, clr_path: str, resolution: int, vmin: float, vmax: float, sep: int, cmap: Any, norm='log', height: int = 5):
        """
        Add Hi-C data to the plot.
        This method reads a Hi-C matrix from a Cooler file, processes it, and adds it to the plot.
        Args:
            clr_path (str): Path to the Cooler file containing Hi-C data.
            resolution (int): Resolution of the Hi-C data in base pairs.
            vmin (float): Minimum value for the color scale.
            vmax (float): Maximum value for the color scale.
            sep (int): Maximum genomic separation (in base pairs) to display on the plot.
            cmap: Colormap to use for visualizing the Hi-C matrix.
            norm: Normalization method to use when displaying the matrix. Defaults to 'log'.
            height (int, optional): Relative height of the Hi-C track. Defaults to 5.
        Notes:
            - The Hi-C matrix is read and padded using `_read_mtx`.
            - The processed matrix is added to the plot using `_add_mtx`.
            - The `trim` parameter is set to `True` when adding the matrix
              to the plot to display it rectangularly instead of triangularly.
        """
        clr = cooler.Cooler(clr_path + f'::/resolutions/{resolution}')
        mtx = self._read_mtx(clr, pad=sep)
        bin_sep = sep // resolution
        self._add_mtx(mtx, vmin, vmax, bin_sep, cmap, norm, height, trim=True)

    def _read_array(self, signal_path: str, binsize: int, mean_subtract=False) -> ArrayLike:
        """
        Reads a signal array from a BigWig file for a specified genomic region.
        Args:
            signal_path (str): Path to the BigWig file containing the signal data.
            binsize (int): Size of the bins to divide the genomic region into.
        Returns:
            numpy.ndarray: Array of signal values for the specified region, with missing values
            filled as NaN. The signal values are summarized using the mean within each bin.
        """

        nbins = (self.end - self.start) // binsize
        with bbi.open(signal_path) as f:
            signal = f.fetch(self.chrom, self.start, self.end, nbins, missing=np.nan, summary='mean')
        if mean_subtract:
            mean = bbi.info(signal_path)['summary'].get('mean', 0)
            signal = signal - mean
        return signal

    def _add_array(self, array: ArrayLike, kind: str, height: int = 2, **kwargs):
        """
        Add an array to the tracks with specified visualization parameters.
        Args:
            array (ArrayLike): The data array to be added as a track.
            kind (str): The type of track to create. Supported values are 'line' and 'fill'.
            height (int, optional): The relative height of the track. Defaults to 2.
            **kwargs: Additional parameters for track customization.
            - For 'line' kind:
                - vmin (float): Minimum value for y axis.
                - vmax (float): Maximum value for y axis.
                - color (str): Color of the line.
            - For 'fill' kind:
                - vmin (float): Minimum value for y axis.
                - vmax (float): Maximum value for y axis.
                - neg_color (str): Color for negative values.
                - pos_color (str): Color for positive values.
        Raises:
            ValueError: If an unsupported track kind is provided.
        """
        if kind == 'line':
            vmin = kwargs['vmin']
            vmax = kwargs['vmax']
            color = kwargs['color']
            track = LineTrack(array, vmin, vmax, color)
        elif kind == 'fill':
            vmin = kwargs['vmin']
            vmax = kwargs['vmax']
            neg_color = kwargs['neg_color']
            pos_color = kwargs['pos_color']
            track = FillTrack(array, vmin, vmax, neg_color, pos_color)
        else:
            raise ValueError(f"Unknown track kind: {kind}")
        self.tracks.append(track)
        self.heights.append(height)
    
    def add_line(self, signal_path: str, binsize: int, vmin: float | None = None, vmax: float | None = None, color: str = 'k', mean_subtract=False, height: int = 2):
        """
        Add a signal as a line track to the region visualization.
        This method reads a signal array from the specified bigwig with a given bin size
        and adds it as a line plot to the visualization.
        Args:
            signal_path (str): Path to the bigwig file to be read.
            binsize (int): Size of the bins to divide the genomic region into.
            vmin (float | None, optional): Minimum value for the y-axis. Defaults to None.
            vmax (float | None, optional): Maximum value for the y-axis. Defaults to None.
            color (str, optional): Color of the line in the plot. Defaults to 'k'.
            mean_subtract (bool, optional): Whether to subtract the genome-wide mean from the signal. Defaults to False.
            height (int, optional): Height of the line plot. Defaults to 2.
        Notes:
            - Missing values are filled as NaN.
            - The signal values are summarized rin bins using the mean.
        """
        signal = self._read_array(signal_path, binsize, mean_subtract=mean_subtract)
        self._add_array(signal, kind='line', height=height, vmin=vmin, vmax=vmax, color=color)
    
    def add_fill(self, signal_path: str, binsize: int, vmin: float | None = None, vmax: float | None = None, neg_color: str = '.8', pos_color: str = 'k', mean_subtract=False, height: int = 2):
        """
        Adds a bigwig signal as a filled track to the region plot.
        This method reads a signal array from the specified bigwig with a given bin size
        and adds it as a filled plot to the visualization.
        Args:
            signal_path (str): Path to the bigwig file.
            binsize (int): Size of the bins to divide the genomic region into.
            vmin (float | None, optional): Minimum value for the y axis. Defaults to None.
            vmax (float | None, optional): Maximum value for the y axis. Defaults to None.
            neg_color (str, optional): Color for negative signal values. Defaults to 'grey'.
            pos_color (str, optional): Color for positive signal values. Defaults to 'k'.
            mean_subtract (bool, optional): Whether to subtract the genome-wide mean from the signal. Defaults to False.
            height (int, optional): Relative height of the filled track. Defaults to 2.
        Notes:
            - Missing values are filled as NaN.
            - The signal values are summarized rin bins using the mean.
        """
        signal = self._read_array(signal_path, binsize, mean_subtract=mean_subtract)
        self._add_array(signal, kind='fill', height=height, vmin=vmin, vmax=vmax, neg_color=neg_color, pos_color=pos_color)
    
    def _read_bed(self, bed_path: str) -> pd.DataFrame:
        """
        Reads a BED-compatible file and returns its contents as a DataFrame.
        Args:
            bed_path (str): The file path to the BED file.
        Returns:
            pandas.DataFrame: A DataFrame containing the contents of the BED file 
            with the schema 'bed6'.
        """
        df = bf.read_table(bed_path, schema='bed6')
        return df

    def add_df(self, df: pd.DataFrame, color: str = 'k', edgecolor: str | None = None, pad: int = 0, height: int = 1):
        """
        Adds a DataFrame with genomic intervals to the region plot with specified color and height.
        Args:
            df (pd.DataFrame): The DataFrame containing genomic data to be added. 
                It is expected to have columns corresponding to chromosome, start, and end positions.
            color (str): The color to be used for the track representation. Defaults to 'k'.
            edgecolor (str, optional): The edge color for the intervals. Defaults to None, will use color if not provided.
            pad (int, optional): The number of base pairs to pad around each interval. Defaults to 0.
            height (int, optional): The relative height of the track in the plot. Defaults to 1.
        """
        subset_df = bf.select(df, (self.chrom, self.start, self.end))
        if edgecolor is None:
            edgecolor = color
        track = IntervalTrack(subset_df, color=color, edgecolor=edgecolor, pad=pad)
        self.tracks.append(track)
        self.heights.append(height)

    def add_df_multi(self, df: pd.DataFrame, hue: str, colors: dict, edgecolors: dict | None = None, pad: int = 0, stratify: bool = False, height: int = 1):
        """
        Adds a DataFrame with genomic intervals with multiple categories to the plot.
        Args:
            df (pd.DataFrame): The input DataFrame containing the data to be plotted.
            hue (str): The column name in the DataFrame used to determine the region category.
            colors (dict): A dictionary mapping categories to their corresponding colors.
            edgecolors (dict, optional): A dictionary mapping categories to their corresponding edge colors.
                Defaults to None, will use colors if not provided.
            pad (int): The number of base pairs to pad around each interval. Defaults to 0.
            stratify (bool): whether to plot each category in its own row. Defaults to False.
            height (int, optional): The relative height of the track to be added. Defaults to 1.
        Raises:
            ValueError: If the specified hue column is not found in the DataFrame.
            ValueError: If colors are not provided for all hue levels present in the DataFrame.
        """
        if hue not in df.columns:
            raise ValueError(f"Hue column '{hue}' not found in DataFrame.")
        missing_hues = set(df[hue].unique()) - set(colors.keys())
        if missing_hues:
            raise ValueError(f"Colors not provided for hue levels: {missing_hues}")
        subset_df = bf.select(df, (self.chrom, self.start, self.end))
        if edgecolors is None:
            edgecolors = colors
        track = MultiIntervalTrack(subset_df, colors=colors, edgecolors=edgecolors, hue=hue, stratify=stratify, pad=pad)
        self.tracks.append(track)
        self.heights.append(height)

    def add_bed(self, bed_path: str, color: str = 'k', edgecolor: str | None = None, pad: int = 0, height: int = 1):
        """
        Adds a BED-compatible file with genomic intervals to the region plot with specified color and height.
        Args:
            bed_path (str): The path to the BED-compatible file with genomic intervals. 
            color (str): The color to be used for the track representation. Defaults to 'k'.
            edgecolor (str, optional): The edge color for the intervals. Defaults to None, will use color if not provided.
            pad (int, optional): The number of base pairs to pad around each interval. Defaults to 0.
            height (int, optional): The relative height of the track in the plot. Defaults to 1.
        """
        df = self._read_bed(bed_path)
        self.add_df(df, color, edgecolor, pad, height)
    
    def add_tads(self, tads: str | pd.DataFrame, facecolor: str = '.6', edgecolor: str = 'k', facealpha: float = .8, edgealpha: Optional[float] = None, orient: Literal["up", "down"] = "up", height: int = 4):
        """
        Adds Topologically Associating Domains (TADs) to the plot.

        This method allows you to overlay TADs on the plot, either by providing
        a file path to a BED file or a pandas DataFrame containing the TAD data.

        Args:
            tads (str | pd.DataFrame): The TAD data to add. Can be a file path to a 
                BED file or a pandas DataFrame with TAD information.
            facecolor (str, optional): The fill color for the TAD regions. Defaults to '.6'.
            edgecolor (str, optional): The color of the edges for the TAD regions. Defaults to 'k'.
            facealpha (float, optional): The transparency level for the TAD fill color. 
                Defaults to 0.8.
            edgealpha (Optional[float], optional): The transparency level for the TAD edge color. 
                If not provided, defaults to the value of `facealpha`.
            orient (Literal["up", "down"], optional): The orientation of the TADs on the plot. 
                Can be "up" or "down". Defaults to "up".
            height (int, optional): The height of the TAD track on the plot. Defaults to 4.

        Raises:
            ValueError: If the `orient` argument is not one of "up" or "down".
        """
        if not isinstance(tads, pd.DataFrame):
            tads = self._read_bed(tads)
        if edgealpha is None:
            edgealpha = facealpha
        subset_df = bf.select(tads, (self.chrom, self.start, self.end))
        track = TADTrack(subset_df, facecolor, edgecolor, facealpha, edgealpha, orient)
        self.tracks.append(track)
        self.heights.append(height)

    def add_ruler(self, loc='bottom', locator: Optional[Any] = None, formatter: Optional[Any] = None, step: Optional[int] = None, fmt_func: Optional[Callable[[int], str]] = None, height: int = 1):
        """
        Adds a ruler track to a plot with customizable location, tick locator, and formatter.
        Args:
            loc (str, optional): The location of the ticks on the ruler. Defaults to 'bottom'.
            locator (Optional[Any], optional): A matplotlib locator instance to determine tick positions. 
                If None, `step` must be provided to create a `MultipleLocator`. Defaults to None.
            formatter (Optional[Any], optional): A matplotlib formatter instance to format tick labels. 
                If None, `fmt_func` must be provided to create a `FuncFormatter`. Defaults to None.
            step (Optional[int], optional): The step size for tick positions if `locator` is not provided. 
                Required if `locator` is None. Defaults to None.
            fmt_func (Optional[Callable[[int], str]], optional): A function to format tick labels if `formatter` 
                is not provided. Required if `formatter` is None. Takes an integer coordinate and returns a string.
                Defaults to None.
        Raises:
            ValueError: If both `locator` and `step` are None.
            ValueError: If both `formatter` and `fmt_func` are None.
        """

        if locator is None:
            if step is None:
                raise ValueError
            locator = plt.MultipleLocator(step)
        if formatter is None:
            if fmt_func is None:
                raise ValueError
            formatter = plt.FuncFormatter(lambda x, pos: fmt_func(x))
        track = RulerTrack(loc, locator, formatter)
        self.tracks.append(track)
        self.heights.append(height)

    def _plot_hic(self, track: MtxTrack, ax: plt.Axes, cax: plt.Axes):
        """
        Plots a Hi-C heatmap on the given axes with specific transformations.
        Args:
            track (MtxTrack): The Hi-C track data to be plotted. Contains attributes such as 
                `data` (the matrix to plot), `vmin` (minimum value for color scaling), 
                `vmax` (maximum value for color scaling), `norm` (normalization for color scaling), 
                `cmap` (colormap), `sep` (separation value for axis limits), and `trim` 
                (boolean indicating whether to trim the plot).
            ax (matplotlib.axes.Axes): The matplotlib axes on which to plot the heatmap.
            cax (matplotlib.axes.Axes): The matplotlib axes for the colorbar.
        """
        sx = 1 / 2 ** .5  # shrink the main diag by sqrt(2)
        sy = 2 ** .5  # increase the height by sqrt(2)
        g = ax.imshow(track.data,
                      interpolation='none',
                      transform=mtransforms.Affine2D().rotate_deg(-45).scale(sx, sy) + ax.transData,
                      vmin=track.vmin,
                      vmax=track.vmax,
                      norm=track.norm,
                      cmap=track.cmap)
        plt.colorbar(g, cax=cax)
        ax.set_ylim(0, track.sep)
        if track.trim:
            ax.set_xlim(track.sep - .5, track.data.shape[0] - track.sep - .5)  # magic numbers to match w/o trimming
        ax.set_aspect(0.5)
        sns.despine()
    
    def _plot_line(self, track: LineTrack, ax: plt.Axes):
        """
        Plots a line track on the given matplotlib axis.
        Args:
            track (LineTrack): The line track object containing data, color, and axis limits.
            ax (matplotlib.axes.Axes): The matplotlib axis on which the line will be plotted.
        """
        ax.plot(track.data, color=track.color)
        ax.set_ylim(track.vmin, track.vmax)
        ax.set_xlim(0, len(track.data))
        sns.despine()
    
    def _plot_fill(self, track: FillTrack, ax: plt.Axes):
        """
        Plots a filled region plot for the given track data on the provided axis.
        This method visualizes the data in the `FillTrack` object by filling the 
        area between the data values and the x-axis. Positive and negative values 
        are filled with different colors.
        Args:
            track (FillTrack): The track object containing the data to be plotted, 
                along with its associated properties such as colors and limits.
            ax (matplotlib.axes.Axes): The matplotlib axis on which the plot will 
                be drawn.
        """
        x = np.arange(len(track.data))
        mask = track.data > 0
        ax.bar(x[mask], track.data[mask], color=track.pos_color, width=1, align='center')
        ax.bar(x[~mask], track.data[~mask], color=track.neg_color, width=1, align='center')
        ax.set_ylim(track.vmin, track.vmax)
        ax.set_xlim(0, len(track.data))
        sns.despine()
    
    def _plot_df(self, track: IntervalTrack, ax: plt.Axes):
        """
        Plots the data from an IntervalTrack object onto the given matplotlib axis.
        If intervals overlap, they are clustered and displayed in separate rows to avoid visual overlap.
        Args:
            track (IntervalTrack): The track containing interval data to be plotted.
                The `track` object must have a `data` attribute (a DataFrame) and a `color`
                attribute specifying the color of the rectangles.
            ax (matplotlib.axes.Axes): The matplotlib axis on which the intervals will be plotted.
        Notes:
            - The y-axis limits are dynamically adjusted based on the number of regions in the largest cluster.
        """
        clusters = bf.cluster(track.data)
        max_regions = 0
        patch_height = .5
        for _, cluster_df in clusters.groupby('cluster'):
            max_regions = max(max_regions, len(cluster_df))
            for i, row in enumerate(cluster_df.itertuples()):
                x, y = row.start - track.pad, -i
                patch_width = row.end - row.start + 2 * track.pad
                patch = Rectangle((x, y), patch_width, patch_height, alpha=.8, edgecolor=track.edgecolor, facecolor=track.color)
                ax.add_patch(patch)
        ax.set_xlim(self.start, self.end)
        ax.set_ylim(-max_regions + patch_height, 2 * patch_height)
    
    def _plot_df_multi(self, track: MultiIntervalTrack, ax: plt.Axes, lax: plt.Axes):
        """
        Plots a MultiIntervalTrack on the provided axes.
        This method visualizes the intervals of a MultiIntervalTrack object by 
        creating rectangular patches for each interval and grouping them by a 
        specified hue. It also generates a legend for the track's color mapping.
        Args:
            track (MultiIntervalTrack): The track object containing interval data, 
                hue information, and color mapping.
            ax (matplotlib.axes.Axes): The main axes on which the intervals are plotted.
            lax (matplotlib.axes.Axes): The legend axes for displaying the color mapping.
        Notes:
            - doesn't deal with overlapping intervals.
        """
        if track.stratify:
            max_regions = track.data[track.hue].nunique()
        else:
            max_regions = 1
        patch_height = .5
        for i, (level, group_df) in enumerate(track.data.groupby(track.hue)):
            if track.stratify:
                plot_level = i
            else:
                plot_level = 0
            color = track.colors[level]
            edgecolor = track.edgecolors[level]
            for _, row in group_df.iterrows():
                x, y = row.start - track.pad, -plot_level
                patch_width = row.end - row.start + 2 * track.pad
                patch = Rectangle((x, y), patch_width, patch_height, alpha=.8, edgecolor=edgecolor, facecolor=color)
                ax.add_patch(patch)
        ax.set_xlim(self.start, self.end)
        ax.set_ylim(-max_regions + patch_height, 2 * patch_height)
        legend_handles = [Patch(facecolor=color, edgecolor=track.edgecolors[level], label=level)
                          for level, color in track.colors.items()]
        lax.legend(handles=legend_handles, loc='center left', borderaxespad=0.)
        lax.set_axis_off()
    
    def _plot_tads(self, track: TADTrack, ax: plt.Axes):
        max_height = 0
        facecolor = mcolors.to_rgba(track.facecolor, track.facealpha)
        edgecolor = mcolors.to_rgba(track.edgecolor, track.edgealpha)
        for tad in track.data.itertuples():
            s = tad.start
            e = tad.end
            if track.orient == 'up':
                height = (e - s) / 2
            elif track.orient == 'down':
                height = - (e - s) / 2
            else:
                raise ValueError(f"Invalid orient value: {track.orient}")
            midpoint = (s + e) / 2
            coords = ((s, 0), (midpoint, height), (e, 0))
            patch = plt.Polygon(coords, facecolor=facecolor, edgecolor=edgecolor)
            ax.add_patch(patch)
            if self.start <= midpoint <= self.end:
                max_height = max(max_height, abs(height))
        ax.set_xlim(self.start, self.end)
        factor = 1.01  # magic number to fit edges
        if track.orient == 'up':
            ax.set_ylim(0, max_height * factor)
        elif track.orient == 'down':
            ax.set_ylim(- max_height * factor, 0)
        ax.set_aspect(1)
        ax.set_axis_off()
    
    def _plot_ruler(self, track: RulerTrack, ax: plt.Axes):
        ax.set_xlim(self.start, self.end)
        ax.xaxis.set_major_locator(track.locator)
        ax.xaxis.set_major_formatter(track.formatter)
        sns.despine()
        ax.yaxis.set_visible(False)
        for key, spine in ax.spines.items():
            spine.set_visible(False)
        if track.loc == 'top':
            ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)

    def plot(self, figsize: tuple[int, int], dpi: int = 100, width_ratios: tuple[int, int] = (40, 1), tick_all: bool = True, tick_step: int = 1_000_000) -> tuple[plt.Figure, ArrayLike]:
        """
        Generate a plot for the tracks with specified layout and formatting.
        Args:
            figsize (tuple[int, int]): The size of the figure in inches (width, height).
            dpi (int): The resolution of the figure in dots per inch. Defaults to 100.
            width_ratios (tuple[int, int], optional): The relative widths of the main plot and the colorbar. Defaults to (40, 1).
            tick_all (bool, optional): Whether to display tick labels for all tracks. Defaults to True.
                If False, will display tick labels only for the bottom track. Ignored if a ruler track is added.
            tick_step (int, optional): The step size for the x-axis ticks in base pairs. Defaults to 1,000,000. Ignored if a ruler track is added.
        Returns:
            tuple: A tuple containing:
                - fig (matplotlib.figure.Figure): The generated figure.
                - axes (numpy.ndarray): The array of axes for the subplots.
        """
        nrows = len(self.tracks)
        ncols = 2
        nticks = (self.end - self.start) // tick_step + 1
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, dpi=dpi, width_ratios=width_ratios, height_ratios=self.heights, squeeze=False)
        ruler_added = False
        for row, track in enumerate(self.tracks):
            ax = axes[row, 0]
            cax = axes[row, 1]
            if isinstance(track, MtxTrack):
                self._plot_hic(track, ax, cax)
            elif isinstance(track, LineTrack):
                self._plot_line(track, ax)
                cax.set_visible(False)
            elif isinstance(track, FillTrack):
                self._plot_fill(track, ax)
                cax.set_visible(False)
            elif isinstance(track, IntervalTrack):
                self._plot_df(track, ax)
                cax.set_visible(False)
            elif isinstance(track, MultiIntervalTrack):
                self._plot_df_multi(track, ax, cax)
            elif isinstance(track, TADTrack):
                self._plot_tads(track, ax)
                cax.set_visible(False)
            elif isinstance(track, RulerTrack):
                self._plot_ruler(track, ax)
                cax.set_visible(False)
                ruler_added = True
            else:
                raise ValueError(f"Unknown track type: {type(track)}")
            
            if not ruler_added:
                ax.xaxis.set_major_locator(plt.LinearLocator(nticks))
                ax.set_xticks(ax.get_xticks())
                if tick_all:
                    ax.set_xticklabels([f"{(self.start + i * tick_step) / 1_000_000:.1f}Mb" for i in range(nticks)])
                else:
                    ax.set_xticklabels([])
        if not ruler_added and not tick_all:
            ax.set_xticks(ax.get_xticks())
            ax.set_xticklabels([f"{(self.start + i * tick_step) / 1_000_000:.1f}Mb" for i in range(nticks)])
        return fig, axes
