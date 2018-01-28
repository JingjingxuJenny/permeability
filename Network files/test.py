from xlrd import open_workbook

from xlutils.copy import copy

rb = open_workbook('pressure_relatveK_4samples.xls',formatting_info=True)
rs = rb.sheet_by_index(0)

wb = copy(rb)

ws = wb.add_sheet('Sheet 2')
ws.write(0, 0, "Appended")

wb.save('example2.xls')
