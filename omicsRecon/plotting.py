import seaborn as sns
import matplotlib.pyplot as plt


def boxplot(draw_df, palette_dict = {},outDir = None, title = None, 
            x=None, y='pair_num', hue=None,figsize = (2.4, 3),
            ylabel=None, dodge=False,legend = False):
    '''
    This function draw boxplot for the input dataframe.
    '''
    if palette_dict == {}:
        palette_dict = ['#EC7063','#3498DB','#63C1A3','#E88AC2','#FB8E68','#9DB0D4']

    if x is None:
        if hue is None:
            x = 'method'
            hue = 'method'
            draw_df['method'] = 'method'
        else:
            x = hue

    plt.figure(figsize=figsize)
    sns.boxplot(x=x, y=y, data=draw_df, hue=hue, dodge=dodge, palette=palette_dict)
    plt.tight_layout()
    plt.xlabel('', size=14)
    plt.ylabel(ylabel, size=14)
    plt.title(title,size = 16)
    if legend:
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    else:
        plt.legend([],[], frameon=False)
    if outDir:
        plt.savefig(f'{outDir}/{title}_boxplot.pdf')