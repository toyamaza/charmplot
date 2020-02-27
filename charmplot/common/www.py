import logging
import os
import subprocess
from datetime import datetime

logger = logging.getLogger(__name__)


def stage_out_plots(folder, vars=[], x=250, y=250):
    www_folder = os.path.join(f"/global/project/projectdirs/atlas/www/{os.environ['USER']}")
    if not os.path.isdir(www_folder):
        os.makedirs(www_folder)
        subprocess.call(['chmod', '755', www_folder])

    # today's date
    today = datetime.now()
    today_string = today.strftime('%Y_%m_%d')
    analysis_folder = os.path.join(www_folder, today_string)
    if not os.path.isdir(analysis_folder):
        os.makedirs(analysis_folder)

    # copy output files over
    subprocess.call(['cp', '-r', folder, analysis_folder])
    subprocess.call(['chmod', '-R', '755', analysis_folder])

    for subdir, dirs, files in os.walk(os.path.join(analysis_folder, folder)):
        for channel in dirs:
            channel_folder = os.path.join(analysis_folder, folder, channel)
            logger.info(f"making index.html for channel {channel}")
            with open(os.path.join(channel_folder, "index.html"), 'w') as index_html:
                index_html.write("<div class=\"ext-box\">\n")
                index_html.write("<TABLE border= \"1\">\n")
                n_files = 0
                # pre-defined variables
                for var in vars:
                    print_line = False
                    for file in os.listdir(channel_folder):
                        if var in file:
                            print_line = True
                            break
                    if not print_line:
                        continue
                    if n_files % 2 == 0:
                        if n_files > 0:
                            index_html.write("</TR>\n")
                        index_html.write("<TR>\n")
                    index_html.write("<TD><CENTER>\n")
                    index_html.write(f"<a href=\"./{channel}_{var}.pdf\">\n")
                    index_html.write(f"<img src=\"{channel}_{var}.png\" height=\"{y}\" width=\"{x}\">\n")
                    index_html.write("</a></CENTER></TD>\n")
                    index_html.write("<TD><CENTER>\n")
                    index_html.write(f"<a href=\"./{channel}_{var}_LOG.pdf\">\n")
                    index_html.write(f"<img src=\"{channel}_{var}_LOG.png\" height=\"{y}\" width=\"{x}\">\n")
                    index_html.write("</a></CENTER></TD>\n")
                    n_files += 1
                # the rest of found variables
                for file in os.listdir(channel_folder):
                    print_line = True
                    for var in vars:
                        if var in file:
                            print_line = False
                    if not print_line:
                        continue
                    if file.endswith(".png") and "_LOG" not in file:
                        if n_files % 2 == 0:
                            if n_files > 0:
                                index_html.write("</TR>\n")
                            index_html.write("<TR>\n")
                        index_html.write("<TD><CENTER>\n")
                        index_html.write(f"<a href=\"./{file.replace('.png', '.pdf')}\">\n")
                        index_html.write(f"<img src=\"{file}\" height=\"{y}\" width=\"{x}\">\n")
                        index_html.write("</a></CENTER></TD>\n")
                        index_html.write("<TD><CENTER>\n")
                        index_html.write(f"<a href=\"./{file.replace('.png', '_LOG.pdf')}\">\n")
                        index_html.write(f"<img src=\"{file.replace('.png', '_LOG.png')}\" height=\"{y}\" width=\"{x}\">\n")
                        index_html.write("</a></CENTER></TD>\n")
                        n_files += 1
