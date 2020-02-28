import logging
import os
import subprocess
from datetime import datetime

logger = logging.getLogger(__name__)


def write_img_line(index_html, file, x, y):
    index_html.write("<TD><CENTER>\n")
    index_html.write(f"<a href=\"./{file}.pdf\">\n")
    index_html.write(f"<img src=\"{file}.png\" height=\"{y}\" width=\"{x}\">\n")
    index_html.write("</a></CENTER></TD>\n")
    index_html.write("<TD><CENTER>\n")
    index_html.write(f"<a href=\"./{file}_LOG.pdf\">\n")
    index_html.write(f"<img src=\"{file}_LOG.png\" height=\"{y}\" width=\"{x}\">\n")
    index_html.write("</a></CENTER></TD>\n")


def place_images(index_html, channel_folder, vars, x, y):
    n_files = 0
    for var in vars:
        if not var.stage_out:
            continue
        print_line = False
        for file in os.listdir(channel_folder):
            if var.name in file:
                print_line = True
                break
        if not print_line:
            continue
        if n_files % 2 == 0 and n_files > 0:
            index_html.write("</TR>\n")
            if var != vars[-1]:
                index_html.write("<TR>\n")
        write_img_line(index_html, f"{os.path.basename(channel_folder)}_{var.name}", x, y)
        n_files += 1


def make_www_folder():
    www_folder = os.path.join(f"/global/project/projectdirs/atlas/www/wcharm")
    if not os.path.isdir(www_folder):
        os.makedirs(www_folder)
        subprocess.call(['chmod', '777', www_folder])
    return www_folder


def make_analysis_folder(www_folder, analysis):
    today = datetime.now()
    today_string = today.strftime('%Y_%m_%d')
    analysis_folder = os.path.join(www_folder, analysis, today_string)
    if not os.path.isdir(analysis_folder):
        os.makedirs(analysis_folder)
    return analysis_folder


def stage_out_plots(analysis, vars=[], x=250, y=250):
    www_folder = make_www_folder()

    analysis_folder = make_analysis_folder(www_folder, analysis)

    # copy output files
    for subdir, dirs, files in os.walk(analysis):
        for channel in dirs:
            subprocess.call(['cp', '-r', os.path.join(analysis, channel), analysis_folder])
    subprocess.call(['chmod', '-R', '755', os.path.join(www_folder, analysis)])

    # create a html.index for each channel
    for subdir, dirs, files in os.walk(analysis_folder):
        for channel in dirs:
            channel_folder = os.path.join(analysis_folder, channel)
            logger.info(f"making index.html for channel {channel}")
            with open(os.path.join(channel_folder, "index.html"), 'w') as index_html:
                index_html.write("<div class=\"ext-box\">\n")
                index_html.write("<TABLE border= \"1\">\n")
                index_html.write("<TR>\n")

                # pre-defined variables
                place_images(index_html, channel_folder, vars, x, y)
