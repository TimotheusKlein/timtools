#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
#
# ad hoc developed custom sqlite3 wrapper for r88_structurama
#
# (C) Timotheus Andreas Klein 2012 - Tim.Klein@gmx.de
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import sqlite3

class sqlctrl():

    def __init__(self, name):
        self.connection = sqlite3.connect(name)
        self.cursor = self.connection.cursor()

    def createindex(self, tname, colname):
        sqlstring = ("CREATE  INDEX IF NOT EXISTS 'main'.'idx_" +
                     tname + "_" +
                     colname + "' ON '" + tname +
                     "' ('" + colname + "' ASC)")
        self.executesql(sqlstring)

    def gettablenames(self):
        sqlstring = ("SELECT name FROM sqlite_master WHERE type = 'table'")
        self.executesql(sqlstring)
        ts = []
        for t in self.cursor.fetchall():
            ts.append(list(t)[0])
        return ts

    def createtable(self, tname, columns, datatypes):
        if tname in self.gettablenames():
            self.droptable(tname)
        sqlstring = ("CREATE TABLE " +
                     tname +
                     " (id INTEGER PRIMARY KEY AUTOINCREMENT,")
        for i, col in enumerate(columns):
            sqlstring = (sqlstring + " '" +
                         str(col) + "' " +
                         str(datatypes[i]) + ",")
        sqlstring = sqlstring.rstrip(",") + ")"
        self.executesql(sqlstring)

    def droptable(self, tname):
        sqlstring = ("DROP TABLE " + tname)
        self.executesql(sqlstring)

    def insertrow(self, tname, columns, data):
        sqlstring = ("INSERT INTO " +
                     tname +
                     " ('")
        for eachcol in columns:
            sqlstring = (sqlstring + str(eachcol) + "', '")
        sqlstring = sqlstring.rstrip(", '") + "') VALUES ("
        for eachentry in data:
            sqlstring += "?, "
        sqlstring = sqlstring.rstrip(", ") + ")"
        try:
            self.cursor.execute(sqlstring, tuple(data))
        except:
            print sqlstring, tuple(data)
            self.cursor.execute(sqlstring, tuple(data))

    def deleterow(self, tname, wheredict):
        sqlstring = ("DELETE FROM " +
                     tname)
        if len(wheredict) > 0:
            sqlstring = sqlstring + " WHERE "
        for k, v in wheredict.iteritems():
            sqlstring = (sqlstring + str(k) + " " +
                             str(v) + " AND ")
        sqlstring = sqlstring.rstrip(" AND ")
        self.executesql(sqlstring)

    def insertcol(self, tname, colname, datatype):
        sqlstring = ("ALTER TABLE " +
                     tname +
                     " ADD COLUMN '" +
                     colname + "' " +
                     datatype)
        self.executesql(sqlstring)

    def changeentry(self, tname, column, newvalue, wheredict):
        sqlstring = ("UPDATE " +
                     tname +
                     " SET '" +
                     column + "' = ? ")
        if len(wheredict) > 0:
            sqlstring = sqlstring + " WHERE "
        for k, v in wheredict.iteritems():
            sqlstring = (sqlstring + str(k) + " " +
                             str(v) + " AND ")
        sqlstring = sqlstring.rstrip(" AND ")
        self.cursor.execute(sqlstring, (newvalue,))

    def selectentries(self, tname, whichcols, wheredict):
        sqlstring = "SELECT "
        if whichcols == '*':
            sqlstring += ' * '
        else:
            for eachcol in whichcols:
                sqlstring = (sqlstring +
                             tname + ".'" + eachcol + "', ")
        sqlstring = (sqlstring.rstrip(", ") +
                     " FROM " + tname)
        if len(wheredict) > 0:
            sqlstring = sqlstring + " WHERE "
            for k, v in wheredict.iteritems():
                sqlstring = (sqlstring + str(k) + " " +
                                 str(v) + " AND ")
            sqlstring = sqlstring.rstrip(" AND ")
        self.executesql(sqlstring)
        t = self.cursor.fetchall()
        t = self.listit(t)
        return t

    def listit(self, t): 
        return list(map(self.listit, t)) if isinstance(t, (list, tuple)) else t 

    def countentries(self, tname, wheredict):
        sqlstring = "SELECT COUNT(*) FROM " + tname
        if len(wheredict) > 0:
            sqlstring = sqlstring + " WHERE "
        for k, v in wheredict.iteritems():
            sqlstring = (sqlstring + str(k) + " " +
                             str(v) + " AND ")
        sqlstring = sqlstring.rstrip(" AND ")
        self.executesql(sqlstring)
        return self.cursor.fetchall()

    def getrow(self, tname, wheredict):
        sqlstring = ("SELECT * FROM " + tname)
        if len(wheredict) > 0:
            sqlstring = sqlstring + " WHERE "
        for k, v in wheredict.iteritems():
            sqlstring = (sqlstring + str(k) + " " +
                             str(v) + " AND ")
        sqlstring = sqlstring.rstrip(" AND ")
        self.executesql(sqlstring)
        return self.cursor.fetchall()

    def getcol(self, tname, colname, wheredict):
        sqlstring = ("SELECT " + tname + ".'" + colname + "' " +
                     " FROM " + tname)
        if len(wheredict) > 0:
            sqlstring = sqlstring + " WHERE "
        for k, v in wheredict.iteritems():
            sqlstring = (sqlstring + str(k) + " " +
                             str(v) + " AND ")
        sqlstring = sqlstring.rstrip(" AND ")
        self.executesql(sqlstring)
        return self.listit(self.cursor.fetchall())

    def getcolnames(self, tname):
        self.cursor.execute("SELECT * FROM " + tname)
        col_name_list = [tuple[0] for tuple in self.cursor.description]
        return col_name_list

    def gettable(self, tname):
        sqlstring = 'SELECT * FROM ' + str(tname)
        self.executesql(sqlstring)
        return self.listit(self.cursor.fetchall())

    def executesql(self, sqlstring):
        trash = self.cursor.execute(sqlstring)

    def finish(self):
        self.connection.commit()
        self.cursor.close()
        
