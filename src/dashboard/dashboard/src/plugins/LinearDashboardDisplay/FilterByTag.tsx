import React, { useState } from 'react'
import { observer } from 'mobx-react'
import { Dialog } from '@headlessui/react'
import {
    Paper, 
    TextField,
    Typography,
    DialogContent,
    DialogActions,
    Button
} from '@mui/material'
import { makeStyles } from 'tss-react/mui'

const useStyles = makeStyles()(theme => ({
    paper: {
      padding: theme.spacing(2),
      margin: theme.spacing(2),
    },
    field: {
      margin: theme.spacing(2),
    },
  }))

const FilterByTagDialog = observer(function (props: {
    model: {
        filterBy: IFilter,
        setFilterBy: (arg: IFilter) => void
    }
    handleClose: () => void
}){
    const { model, handleClose } = props
    const { classes } = useStyles()
    const { filterBy } = model
    const [flagInclude, setFlagInclude] = useState(filterBy.flagInclude)
    const [flagExclude, setFlagExclude] = useState(filterBy.flagExclude)
    const [tag, setTag] = useState(filterBy.tagFilter?.tag || '')
    const [tagValue, setTagValue] = useState(filterBy.tagFilter?.value || '')
    const [readName, setReadName] = useState(filterBy.readName || '')
    const validTag = tag.match(/^[A-Za-z][A-Za-z0-9]$/)

    return (
    <Dialog open onClose={handleClose} title="Filter options">
    <DialogContent>

        <Paper className={classes.paper} variant="outlined">
          <Typography>
            Filter by tag name and value. Use * in the value field to get all
            reads containing any value for that tag. Example: filter tag name SA
            with value * to get all split/supplementary reads. Other examples
            include HP for haplotype, or RG for read group
          </Typography>
          <TextField
            className={classes.field}
            value={tag}
            onChange={event => setTag(event.target.value)}
            placeholder="Enter tag name"
            inputProps={{ maxLength: 2 }}
            error={tag.length === 2 && !validTag}
            helperText={tag.length === 2 && !validTag ? 'Not a valid tag' : ''}
          />
          <TextField
            className={classes.field}
            value={tagValue}
            onChange={event => setTagValue(event.target.value)}
            placeholder="Enter tag value"
          />
        </Paper>
        <DialogActions>
          <Button
            variant="contained"
            color="primary"
            autoFocus
            type="submit"
            onClick={() => {
              model.setFilterBy({
                flagInclude,
                flagExclude,
                readName,
                tagFilter:
                  tag !== ''
                    ? {
                        tag,
                        value: tagValue,
                      }
                    : undefined,
              })
              handleClose()
            }}
          >
            Submit
          </Button>
          <Button
            variant="contained"
            color="secondary"
            onClick={() => handleClose()}
          >
            Cancel
          </Button>
        </DialogActions>
      </DialogContent>
    </Dialog>
)})

export interface IFilter {
    flagExclude: number
    flagInclude: number
    readName?: string
    tagFilter?: { tag: string; value: string }
    filterReads?: any
  }

export default FilterByTagDialog