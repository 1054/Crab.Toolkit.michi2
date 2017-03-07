;
; Updated: 2013-03-27
; 
FUNCTION CrabTableReadItem, FilePath, ColumnHeader, RowText

    ThisColumn = CrabTableReadColumn(FilePath, ColumnHeader)
    ThisRow = CrabTableReadRow(FilePath, RowText, RowID=ThisRowID)
    
    RETURN, ThisColumn[ThisRowID]
    
END