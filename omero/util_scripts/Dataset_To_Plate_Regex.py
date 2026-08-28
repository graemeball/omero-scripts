#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
#   Copyright (C) 2006-2026 University of Dundee. All rights reserved.
#
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program; if not, write to the Free Software Foundation, Inc.,
#   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# ------------------------------------------------------------------------------


"""
Adds a Dataset of Images to a new Plate, extracting Well info from Image names.
E.g.
1. Well Row and Column only:
   - image_name='WellB01_image.tif'
   - regex='Well(?P<row>[A-P])(?P<col>[0-9]{2})_image.tif'
2. Well Row and Column plus Field number (i.e. position/site/series):
   - image_name='WellB01_WT_Pos01.tif'
   - regex='Well(?P<row>[A-P])(?P<col>[0-9]{2})_WT_Pos(?P<field>[0-9]{2}).tif'

Optionally adds the new Plate to a new or existing Screen.

Notes:
- script assumes <row> is alphabetical (A-P), <col> is integer (1-24)

See http://help.openmicroscopy.org/scripts.html
"""

# @authors Graeme Ball, Will Moore
# <a href="mailto:g.ball@dundee.ac.uk">g.ball@dundee.ac.uk</a>
# @version 1.0

# Use & misuse scenarios tested (v1.0):
# 1a. multiple datasets specified - Pass (processed first)
# 1b. dataset id not found -  Pass (message: Dataset not found)
# 2a. regex invalid / cannot compile - Pass (message: regex error)
# 2b. regex doesn't match - Pass (message: 0 added, stdout: no match)
# 2c. invalid row, col matched by regex - Pass
#      (skips invalid, stdout: rol and col indices not in valid lists)
# 3a. handle <field> pattern not an integer - Pass (stdout: invalid)
# 3b. multiple images per well, field pattern present in some - Pass
# 3c. multiple images per well, no field pattern - Pass
# 3d. multiple images per well, field pattern present in all - Pass
# 3e. handle multiple fields with same number - Pass (all added)
# 4a. name a new screen - Pass
# 4b. add to existing screen by id (exists) - Pass
# 4c. add to existing screen by id (doesn't exist) - Pass
#     (Plate created, message: no such Screen)
# 5a. option to not remove images from dataset - Pass
# 5b. and then try to add images to a plate again - Pass
#      (refuses until links removed, links shown in stdout)

import omero.scripts as scripts
from omero.gateway import BlitzGateway
import omero

from omero.rtypes import rint, rlong, rstring, robject

import re
import string


def extract_well_row_col_field(image_name, regex_compiled):
    """
    Return tuple of ((row, column, field), info)
      where: (row, column, field)
      is 0-based index tuple extracted from image_name and Field can be None
      using regex with named <row>, <col> and (optional) <field> patterns.
    E.g. image_name='WellB01_image.tif',
      pattern='Well(?P<row>[A-P])(?P<col>[0-1]{1,2})_image.tif'
    N.B. assumes <row> is alphabetical (A-P), <col> is integer (1-24);
    Returns (None, info) if pattern does not match or row/col not found.
    """
    # row=A-P,col=1-24; i.e. maximum 384 well plate
    row_labels = list(string.ascii_uppercase)[0:16]
    col_labels = list(range(1, 25, 1))

    info = ""

    # search for regex pattern, convert to 0-based row,col,field indices
    m = regex_compiled.search(image_name)
    info += f"{image_name}: "
    if m is None:
        info += "No match! "
        return None, info
    else:
        mg = m.groupdict()
        print(f"m.groupdict()={m.groupdict()}")  # FIXME, remove
        try:
            row = mg['row'].upper()
            col = int(mg['col'])
        except KeyError:
            info += "<row> and/or <col> named matches missing. "
            return None, info
        if 'field' in mg.keys():
            field = mg['field']
            try:
                field = int(field)
            except Exception as e:
                info += f"<field> should match a number - {e}"
                field = None
        else:
            field = None

        try:
            row_index = row_labels.index(row)
        except ValueError:
            info += f"row index not in {row_labels}. "
            row_index = None

        try:
            col_index = col_labels.index(col)
        except ValueError:
            info += f"column index not in {col_labels}. "
            col_index = None

        if (row_index is not None) and (col_index is not None):
            return (row_index, col_index, field), info
        else:
            return None, info


def add_fields_to_plate_well(conn, image_field_list, plate_id,
                             row, column, remove_from=None):
    """
    Add list of (Image, Field) to a Plate, creating a new Well
    at the specified row and column.
    The order of Image addition is according to Field if defined
    (otheriwse None goes last).
    NB - This will fail if the Well already exists.
    """
    update_service = conn.getUpdateService()

    well = omero.model.WellI()
    well.plate = omero.model.PlateI(plate_id, False)
    well.column = rint(column)
    well.row = rint(row)

    def sort_list_of_ab_by_b(lab):
        # sort list of (a, b) according to b, where b=None goes last
        return sorted(lab, key=lambda lab: (lab[1] is None, lab[1]))

    image_field_list_sorted = sort_list_of_ab_by_b(image_field_list)

    images = []
    try:
        for image, field in image_field_list_sorted:
            ws = omero.model.WellSampleI()
            ws.image = omero.model.ImageI(image.id, False)
            ws.well = well
            well.addWellSample(ws)
            images.append(image)
        update_service.saveObject(well)
    except Exception as e:
        print(f"Exception adding Images to row,col={row},{column}: {e}")
        return 0

    # remove from Dataset
    for image in images:
        if remove_from is not None:
            links = list(image.getParentLinks(remove_from.id))
            link_ids = [link.id for link in links]
            conn.deleteObjects('DatasetImageLink', link_ids)

    return len(image_field_list)


def dataset_to_plate(conn, script_params):
    """Try to add Dataset Images to Plate using parameters provided."""

    message = ""

    # get script parameters and update service
    dtype = script_params['Data_Type']
    dataset_ids = script_params['IDs']
    dataset_id = dataset_ids[0]
    if len(dataset_ids) > 1:
        print(f"Processing first Dataset only (id={dataset_id})")
    well_info_regex = script_params['Well_Info_Regex']
    remove_from_dataset = script_params['Remove_From_Dataset']
    screen_id = None
    if "Screen" in script_params:
        screen = script_params["Screen"]
        try:
            screen_id = int(screen)
        except ValueError:
            pass  # i.e. screen is a string and screen_id is None
    else:
        screen = None

    update_service = conn.getUpdateService()

    # check the regex compiles without error!
    try:
        regex_compiled = re.compile(well_info_regex)  # raises RegexError
    except re.error as e:
        message += f"Error! for regex '{well_info_regex}': {e}"
        return None, message

    # Get dataset using ID - abort if Wells already linked / no permission
    dataset = conn.getObject(dtype, dataset_id)
    if dataset is None:
        message += f"Dataset {dataset_id} not found! "
        return None, message

    def has_images_linked_to_well(dataset):
        params = omero.sys.ParametersI()
        query = "select img, well from Well as well "\
                "left outer join well.wellSamples as ws " \
                "left outer join ws.image as img "\
                "where img.id in (:ids)"
        params.addIds([i.getId() for i in dataset.listChildren()])
        imgs_wells = conn.getQueryService().projection(
            query, params, conn.SERVICE_OPTS)
        if len(imgs_wells) > 0:
            for img, well in imgs_wells:
                image = img.getValue()
                print(f"Image:{image.id.val} ({image.name.val}) is linked " +
                      "to Well {well.getValue().id.val}")
            return True

    if has_images_linked_to_well(dataset):
        print(f"Dataset {dataset_id} already has image-well links! ")
        return None, message

    if not dataset.canLink():
        print(f"No permission to add images from dataset {dataset_id}. ")
        return None, message

    # Do we try to remove images from Dataset and Delte Datset when/if empty?
    remove_from = None
    if remove_from_dataset:
        remove_from = dataset

    # build a dict of Image lists for Wells identified in Dataset using regex
    # dict keys are (row, col) tuples, values are Lists of (Image, field)
    images = list(dataset.listChildren())
    plate_well_images = {}
    for image in images:
        row_col_field, info = extract_well_row_col_field(
            image.name, regex_compiled)
        print(info)
        if row_col_field is not None:
            print(f"row,col,field={row_col_field}. ")
            row_col = row_col_field[0:2]
            field = row_col_field[2]
            if row_col not in plate_well_images:
                plate_well_images[row_col] = []
            plate_well_images[row_col].append((image, field))

    # find Screen if specified by ID or create new Screen if Name provided
    newscreen = None
    if screen_id:
        screen = conn.getObject("Screen", screen_id)
        if screen is None:
            message += f"No such Screen (id={screen_id}). "
    elif screen is not None:
        newscreen = omero.model.ScreenI()
        newscreen.name = rstring(screen)
        newscreen = update_service.saveAndReturnObject(newscreen)
        screen = conn.getObject("Screen", newscreen.getId().getValue())
        screen_id = screen.id
        screen_name = screen.name
        print(f"Created new Screen '{screen_name}' (id={screen_id}). ")

    # create Plate & link to Screen if specified
    plate = omero.model.PlateI()
    plate.name = omero.rtypes.RStringI(dataset.name)
    plate.columnNamingConvention = rstring(str('number'))
    plate.rowNamingConvention = rstring(str('letter'))
    plate = update_service.saveAndReturnObject(plate)
    plate_name = plate.getName().getValue()
    plate_id = plate.getId().getValue()
    print(f"New Plate created: {plate_name} (id={plate_id}). ")
    if screen is not None:
        if screen.canLink():
            link = omero.model.ScreenPlateLinkI()
            link.parent = omero.model.ScreenI(screen.id, False)
            link.child = omero.model.PlateI(plate.id.val, False)
            update_service.saveObject(link)
            message += f"Linked new Plate to Screen (id={screen_id}). "
        else:
            message += "Could not link Plate to Screen! "
            message += f"(screen.id={screen_id}). "

    # iterate over image_field_list adding Images to Plate one Well at a time
    added_count = 0
    for row_col, image_field_list in plate_well_images.items():
        added_count += add_fields_to_plate_well(
            conn, image_field_list, plate.getId().getValue(),
            row_col[0], row_col[1], remove_from)

    message += f"Added {added_count} Images to Plate (id={plate_id}). "

    if newscreen is not None:
        robj = newscreen
    elif plate is not None:
        robj = plate
    else:
        robj = None

    return robj, message


def run_script():
    """
    The main entry point of the script, as called by the client via the
    scripting service, passing the required parameters.
    """

    data_types = [rstring('Dataset')]

    client = scripts.client(
        "Dataset_To_Plate_Regex.py",

        """Take all Images found in a Dataset and add them to a new Plate,
        extracting Well row, column and optionally field info
        from Image names using a python regex (regular expression).
        E.g. for Images 'WellB1_WT_Pos002.tif', 'WellC02_WT_series001.tif'
        a matching regex would be:
        'Well(?P<row>[A-P])(?P<col>[0-9]{1,2})_WT_Pos(?P<field>[0-9]{3}).tif'
        N.B. named <row> and <col> patterns must be captured by the regex!
        --
        Optionally add the Plate to a new or existing Screen.
        For help see:
        - OMERO scripts: http://help.openmicroscopy.org/scripts.html
        - Regular expressions:
        https://cellprofiler-manual.s3.amazonaws.com/CPmanual/Metadata.html""",

        scripts.String(
            "Data_Type", optional=False, grouping="1",
            description="Choose source of images (only Dataset supported)",
            values=data_types, default="Dataset"),

        scripts.List(
            "IDs", optional=False, grouping="2",
            description="Dataset ID to convert to new Plate (N.B. one only)."
        ).ofType(rlong(0)),

        scripts.String(
            "Well_Info_Regex", optional=False, grouping="3", default="",
            description="Regex to capture <row>, <col> and optionally " +
            "<field> from image names."),

        scripts.String(
            "Screen", grouping="4",
            description="Option: put Plate in a Screen. " +
            "Enter ID of existing screen " +
            "or Name of new Screen"),

        scripts.Bool(
            "Remove_From_Dataset", grouping="5", default=True,
            description="Remove Images from Dataset as added to Plate"),

        version="1.0",
        authors=["Graeme Ball", "William Moore", "OME Team"],
        institutions=["University of Dundee"],
        contact="g.ball@dundee.ac.uk",
    )

    try:
        script_params = client.getInputs(unwrap=True)

        # wrap client to use the Blitz Gateway
        conn = BlitzGateway(client_obj=client)

        # Convert Dataset to Plate. Returns new plate or screen.
        new_obj, message = dataset_to_plate(conn, script_params)

        client.setOutput("Message", rstring(message))
        if new_obj is not None:
            client.setOutput("New_Object", robject(new_obj))

    finally:
        client.closeSession()


if __name__ == "__main__":
    run_script()
